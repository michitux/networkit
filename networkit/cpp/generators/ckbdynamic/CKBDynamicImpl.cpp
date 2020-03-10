#include "CKBDynamicImpl.h"
#include "CommunityDeathEvent.h"
#include "CommunityBirthEvent.h"
#include "CommunitySplitEvent.h"
#include "CommunityMergeEvent.h"
#include "CustomCommunitySizeDistribution.h"
#include "CustomCommunityMembershipDistribution.h"
#include "PowerlawCommunitySizeDistribution.h"
#include "PowerlawCommunityMembershipDistribution.h"
#include "../../auxiliary/SignalHandling.h"
#include "../../auxiliary/Timer.h"
#include <tlx/unused.hpp>

namespace NetworKit {
	namespace CKBDynamicImpl {
		void CKBDynamicImpl::addEdge(node u, node v, bool nodeJoined) {
			auto e = Community::canonicalEdge(u, v);
			index ts = currentTimeStep;

			if (edgeSharpness < 1 && nodeJoined && currentTimeStep > 0) {
				count offset = edge_sharpness_distribution(urng);
				if (offset < ts) {
					ts -= offset;
				} else {
					ts = 0;
				}
			}

			eventStream.addEdge(ts, e.first, e.second);
		}

		void CKBDynamicImpl::removeEdge(node u, node v, bool nodeLeft) {
			auto e = Community::canonicalEdge(u, v);
			index ts = currentTimeStep;

			if (edgeSharpness < 1 && nodeLeft && currentTimeStep > 0) {
				count offset = edge_sharpness_distribution(urng);
				if (offset + ts < numTimesteps) {
					ts += offset;
				} else {
					ts = numTimesteps;
				}
			}

			eventStream.removeEdge(ts, e.first, e.second);
		}

		void CKBDynamicImpl::addNodeToCommunity(node u, CommunityPtr com) {
			assert(hasNode(u));
			if (com != globalCommunity) {
				if (desiredMemberships[u] == nodeCommunities[u].size()) {
					nodesWithOverassignments.insert(u);
				}
				nodeCommunities[u].insert(com);
				eventStream.nodeJoinsCommunity(currentTimeStep, u, com->getId());
				++currentCommunityMemberships;
				if (desiredMemberships[u] == nodeCommunities[u].size()) {
					nodesWantingAssignments.erase(u);
				}
				if (com->getNumberOfNodes() == com->getDesiredNumberOfNodes()) {
					communitiesWantingMembers.erase(com);
				}

				assert(communitiesWantingMembers.contains(com) == (com->getNumberOfNodes() < com->getDesiredNumberOfNodes()));
			}
		}

		void CKBDynamicImpl::removeNodeFromCommunity(node u, CommunityPtr com) {
			assert(hasNode(u));
			if (com != globalCommunity) {
				if (desiredMemberships[u] == nodeCommunities[u].size()) {
					nodesWantingAssignments.insert(u);
				}
				// The node has already been removed from the community
				if (com->getNumberOfNodes() + 1 == com->getDesiredNumberOfNodes()) {
					communitiesWantingMembers.insert(com);
				}
				nodeCommunities[u].erase(com);
				if (desiredMemberships[u] == nodeCommunities[u].size()) {
					nodesWithOverassignments.erase(u);
				}
				eventStream.nodeLeavesCommunity(currentTimeStep, u, com->getId());
				--currentCommunityMemberships;
				assert(communitiesWantingMembers.contains(com) == (com->getNumberOfNodes() < com->getDesiredNumberOfNodes()));
			}
		}

		void CKBDynamicImpl::addCommunity(CommunityPtr com) {
			if (com->isAvailable()) {
				availableCommunities.insert(com);
			} else {
				availableCommunities.erase(com);
			}
			communities.insert(com);
			// The global community is not an issue here as it will always have 0 desired members
			if (com->getNumberOfNodes() < com->getDesiredNumberOfNodes()) {
				communitiesWantingMembers.insert(com);
			}
		}

		void CKBDynamicImpl::removeCommunity(CommunityPtr com) {
			assert(com->getNumberOfNodes() == 0);
			availableCommunities.erase(com);
			communities.erase(com);
			communitiesWantingMembers.erase(com);
		}

		void CKBDynamicImpl::desiredSizeChanged(CommunityPtr com) {
			assert(com != globalCommunity);

			if (com->getNumberOfNodes() < com->getDesiredNumberOfNodes()) {
				communitiesWantingMembers.insert(com);
			} else {
				communitiesWantingMembers.erase(com);
			}
		}

		index CKBDynamicImpl::nextCommunityId() {
			index result = maxCommunityId;
			++maxCommunityId;
			return result;
		}

		index CKBDynamicImpl::drawIndex(index a, index b) {
			std::uniform_int_distribution<index>::param_type p(a, b - 1);
			return uniform_distribution(urng, p);
		}

		index CKBDynamicImpl::drawIndex(index b) {
			return drawIndex(0, b);
		}

		double CKBDynamicImpl::drawBinomial(count numTrials, double probability) {
			std::binomial_distribution<count>::param_type p(numTrials, probability);
			return binomial_distribution(urng, p);
		}

		double CKBDynamicImpl::drawProbability() {
			return random_probability_distribution(urng);
		}

		CKBDynamicImpl::CKBDynamicImpl(const CKBDynamic::param_type &params) :
			urng(Aux::Random::integer()),
			communitySizeSampler(nullptr),
			membershipDistribution(nullptr),
			random_probability_distribution(0, 1),
			edge_sharpness_distribution(params.edgeSharpness),
			maxCommunityId(0),
			sumOfDesiredMemberships(0),
			currentTimeStep(0),
			eventStream(params.numTimesteps + 1),
			n(params.n),
			communityEventProbability(params.communityEventProbability),
			nodeEventProbability(params.nodeEventProbability),
			perturbationProbability(params.perturbationProbability),
			epsilon(params.epsilon),
			edgeSharpness(params.edgeSharpness),
			tEffect(params.tEffect),
			numTimesteps(params.numTimesteps),
			currentCommunityMemberships(0) {

			if (params.G != nullptr && params.C != nullptr) {
				communitySizeSampler.reset(new CustomCommunitySizeDistribution(*params.G, *params.C));
				epsilon = static_cast<CustomCommunitySizeDistribution*>(communitySizeSampler.get())->getEpsilon();
				membershipDistribution.reset(new CustomCommunityMembershipDistribution(*params.G, *params.C));

			} else {
				communitySizeSampler.reset(new PowerlawCommunitySizeDistribution(params.minCommunitySize, params.maxCommunitySize, params.communitySizeExponent, params.intraCommunityEdgeProbability, params.intraCommunityEdgeExponent));
				membershipDistribution.reset(new PowerlawCommunityMembershipDistribution(params.minCommunityMembership, params.maxCommunityMembership, params.communityMembershipExponent));
			}

			double expectedNumberOfCommunities = membershipDistribution->getAverageMemberships() * n / communitySizeSampler->getAverageSize();
			if (expectedNumberOfCommunities < membershipDistribution->getMaximumMemberships()) {
				throw std::runtime_error("Error: Graph impossible to realize, in expectation, there will be " + std::to_string(expectedNumberOfCommunities) + " communities but there may be a node that wants to be part of " + std::to_string(membershipDistribution->getMaximumMemberships()) + " communities.");
			}

		}

		std::vector<GraphEvent> CKBDynamicImpl::getGraphEvents() {
			this->assureFinished();
			return eventStream.getGraphEvents();
		}

		std::vector<CommunityEvent> CKBDynamicImpl::getCommunityEvents() {
			this->assureFinished();
			return eventStream.getCommunityEvents();
		}

		void CKBDynamicImpl::generateNode() {
			node u = desiredMemberships.size();
			count m = membershipDistribution->drawMemberships();
			desiredMemberships.push_back(m);
			for (count i = 0; i < m; ++i) {
				nodeMemberships.push_back(u);
			}
			sumOfDesiredMemberships += m;
			nodesAlive.insert(u);
			nodeCommunities.emplace_back();
			globalCommunity->addNode(u);
			if (m > 0) {
				nodesWantingAssignments.insert(u);
			}
			eventStream.addNode(currentTimeStep, u);
		}

		void CKBDynamicImpl::eraseNode() {
			node u = nodesAlive.at(drawIndex(nodesAlive.size()));
			sumOfDesiredMemberships -= desiredMemberships[u];
			desiredMemberships[u] = 0;

			while (nodeCommunities[u].size() > 0) {
				CommunityPtr com = *nodeCommunities[u].begin();
				com->removeNode(u);
			}

			assert(nodesAlive.contains(u));
			globalCommunity->removeNode(u);
			nodesWantingAssignments.erase(u);
			nodesAlive.erase(u);
			eventStream.removeNode(currentTimeStep, u);
		}

		void CKBDynamicImpl::run() {
			if (hasRun) throw std::runtime_error("Error, run has already been called");

			Aux::SignalHandler handler;

			// initialization
			globalCommunity = CommunityPtr(new Community(*this));
			globalCommunity->changeEdgeProbability(epsilon);
			communities.erase(globalCommunity);
			availableCommunities.erase(globalCommunity);
			currentTimeStep = 0;

			for (node u = 0; u < n; ++u) {
				generateNode();
			}

			const count initialNumberOfNodes = nodesAlive.size();

			count sumOfDesiredMembers = 0;

			while (sumOfDesiredMembers < sumOfDesiredMemberships) {
				handler.assureRunning();
				count communitySize = communitySizeSampler->drawCommunitySize();;

				CommunityPtr com(new Community(*this));
				com->setDesiredNumberOfNodes(communitySize);
				sumOfDesiredMembers += communitySize;
			}

			assignNodesToCommunities();

			double deathProbability = 0.25, birthProbability = 0.25, splitProbability = 0.25, mergeProbability = 0.25;
			tlx::unused(mergeProbability);

			for (currentTimeStep = 1; currentTimeStep <= numTimesteps; ++currentTimeStep) {
				handler.assureRunning();
				const count numCommunityEvents = drawBinomial(communities.size(), communityEventProbability);

				const count numNodeEvents = drawBinomial(communities.size(), nodeEventProbability);

				INFO("Timestep ", currentTimeStep, " generating ", numCommunityEvents, " community events and ", numNodeEvents, " node events");

				for (count i = 0; i < numCommunityEvents; ++i) {
					{ // adjust event probabilities
						const double x = sumOfDesiredMemberships * 1.0 / sumOfDesiredMembers;
						splitProbability = birthProbability = 0.5 * x / (1 + x);
						mergeProbability = deathProbability = 0.5 - birthProbability;
					}
					handler.assureRunning();
					double r = drawProbability();
					if (r < birthProbability) {
						// generate new community
						count coreSize = communitySizeSampler->getMinSize();
						count targetSize = communitySizeSampler->drawCommunitySize();
						sumOfDesiredMembers += targetSize;
						currentEvents.emplace_back(new CommunityBirthEvent(coreSize, targetSize, tEffect, *this));
					} else if (r < birthProbability + deathProbability) {
						// let a community die
						if (availableCommunities.size() > 0) {
							CommunityPtr com = availableCommunities.at(drawIndex(availableCommunities.size()));
							sumOfDesiredMembers -= com->getDesiredNumberOfNodes();
							count coreSize = communitySizeSampler->getMinSize();
							currentEvents.emplace_back(new CommunityDeathEvent(com, coreSize, tEffect, *this));
							assert(!com->isAvailable());
						} else {
							WARN("No community available for death event.");
						}
					} else if (r < birthProbability + deathProbability + splitProbability) {
						// Split a community
						if (availableCommunities.size() > 0) {
							CommunityPtr com = availableCommunities.at(drawIndex(availableCommunities.size()));
							sumOfDesiredMembers -= com->getDesiredNumberOfNodes();
							count comSizeA = communitySizeSampler->drawCommunitySize();
							sumOfDesiredMembers += comSizeA;
							count comSizeB = communitySizeSampler->drawCommunitySize();
							sumOfDesiredMembers += comSizeB;
							currentEvents.emplace_back(new CommunitySplitEvent(com, comSizeA, comSizeB, tEffect, *this));
							assert(!com->isAvailable());
						} else {
							WARN("No community available for splitting.");
						}
					} else {
						// merge two communities
						if (availableCommunities.size() > 1) {
							index ia = drawIndex(availableCommunities.size());
							index ib = drawIndex(1, availableCommunities.size());
							if (ia == ib) {
								ib = 0;
							}

							CommunityPtr comA = availableCommunities.at(ia);
							sumOfDesiredMembers -= comA->getDesiredNumberOfNodes();
							CommunityPtr comB = availableCommunities.at(ib);
							sumOfDesiredMembers -= comB->getDesiredNumberOfNodes();

							count targetSize = communitySizeSampler->drawCommunitySize();
							sumOfDesiredMembers += targetSize;
							currentEvents.emplace_back(new CommunityMergeEvent(comA, comB, targetSize, tEffect, *this));
							assert(!comA->isAvailable());
							assert(!comB->isAvailable());
						} else {
							WARN("No two communities available for merge.");
						}
					}
				} // generated all new community events

				// generate node events
				const double wantedNodeFraction = initialNumberOfNodes * 1.0 / nodesAlive.size();
				const double nodeBirthProbability = wantedNodeFraction / (1 + wantedNodeFraction);

				// First generate all death events, then all birth events.
				// This ensures that no node that is born in this time step dies again in this time step.
				const count nodesBorn = drawBinomial(numNodeEvents, nodeBirthProbability);

				for (count j = 0; j < (numNodeEvents - nodesBorn); ++j) {
					eraseNode();
				}
				for (count j = 0; j < nodesBorn; ++j) {
					generateNode();
				}

				// Trigger all current events
				for (size_t e = 0; e < currentEvents.size();) {
					handler.assureRunning();
					currentEvents[e]->nextStep();

					if (!currentEvents[e]->isActive()) {
						std::swap(currentEvents[e], currentEvents.back());
						currentEvents.pop_back();
					} else {
						++e;
					}
				}

				// edge perturbations
				if (perturbationProbability > 0) {
					globalCommunity->perturbEdges(perturbationProbability);

					for (CommunityPtr com : communities) {
						handler.assureRunning();
						com->perturbEdges(perturbationProbability);
					}
				}

				assignNodesToCommunities();

				// adjust event probabilities
				{
					INFO("Current memberships: ", currentCommunityMemberships, " desired: ", sumOfDesiredMemberships, ", desired members after events: ", sumOfDesiredMembers," number of communities: ", communities.size(), " available: ", availableCommunities.size(), " active events ", currentEvents.size());
				}
			}

			availableCommunities.clear();
			communities.clear();
			nodeCommunities.clear();
			globalCommunity = nullptr;
			currentEvents.clear();

			eventStream.run();

			hasRun = true;
		}

		void CKBDynamicImpl::assignNodesToCommunities() {
			count totalMissingMembers = 0;


			// Step 1: determine set of nodes wanting communities and communities wanting nodes

			Aux::Timer timer;
			timer.start();

			for (const CommunityPtr& com : communitiesWantingMembers) {
				const count desired = com->getDesiredNumberOfNodes();
				assert(desired >= communitySizeSampler->getMinSize());
				const count actual = com->getNumberOfNodes();
				assert(actual < desired);

				totalMissingMembers += desired - actual;
			}

			if (totalMissingMembers == 0) return;

			count totalMissingMemberships = 0;

			struct NodeMembership {
				count desired;
				count toAssign;
			};

			robin_hood::unordered_flat_map<node, NodeMembership> nodesWithMissingCommunities;

			for (node u : nodesWantingAssignments) {
				const count desired = desiredMemberships[u];
				const count actual = nodeCommunities[u].size();

				assert(nodesAlive.contains(u));
				assert(desired > actual);

				count toAssign = desired - actual;
				totalMissingMemberships += toAssign;
				nodesWithMissingCommunities.insert({u, NodeMembership{desired, toAssign}});
			}

#ifndef NDEBUG
			for (node u : nodesAlive) {
				const count desired = desiredMemberships[u];
				const count actual = nodeCommunities[u].size();

				assert(nodesWithMissingCommunities.contains(u) == (desired > actual));
			}
#endif

			timer.stop();
			INFO("Needed ", timer.elapsedMicroseconds() * 1.0 / 1000, "ms to collect initial candidates, ", totalMissingMembers, " members to be found, ", totalMissingMemberships, " memberships wanted");

			// Step 2: under-assignment: if nodes want more members than communities want members
			// Step 2a: try removing nodes from communities to compensate mismatch, i.e., let communities want more members
			// If we have more node slots free than
			// communities that are missing some members,
			// try to find nodes that got additional
			// members where we can remove some of them.
			// In order to not to disturb merge and split
			// events, we do not take communities that are
			// not available.
			// FIXME: make split and merge events more
			// robust for random node exchanges.
			if (totalMissingMembers < totalMissingMemberships) {
				timer.start();
				for (size_t i = 0; i < nodesWithOverassignments.size() && totalMissingMembers < totalMissingMemberships;) {
					node u = nodesWithOverassignments.sample_item(i);

					assert(nodeCommunities[u].size() > desiredMemberships[u]);

					bool reducedToDesired = false;
					for (size_t ci = 0; ci < nodeCommunities[u].size() && totalMissingMembers < totalMissingMemberships;) {
						const CommunityPtr &com = nodeCommunities[u].sample_item(ci);
						if (com->canRemoveNode()) {
							com->removeNode(u);
							++totalMissingMembers;

							if (nodeCommunities[u].size() == desiredMemberships[u]) {
								reducedToDesired = true;
								break;
							}
						} else {
							// only increment if we did not remove this community
							// if we removed this community, the item at position ci has been replaced
							// and we can sample again at position ci
							++ci;
						}
					}

					if (!reducedToDesired) {
						// only increment if we did not remove all overassignments
						// if we removed all overassignments, the node has been removed from nodesWithOverAssignments
						// and we can sample again at position i
						++i;
					}
				}
				timer.stop();
				INFO("Needed ", timer.elapsedMicroseconds() * 1.0 / 1000, "ms to remove additional nodes from communities, now wanting ", totalMissingMembers, " members");
			}

			// Step 2b: if this did not succeed:
			/* Get a vector with all nodes wanting assignments in the respective multiplicity.
			   Sample from this vector totalMissingMembers nodes using fisher-yates-shuffle.
			   However, first ensure that every node has at least one member and otherwise include them in the sample.
			   If this sample is still too big -> ? - orphans? let nodes die? force-remove some nodes with more than one member?
			*/
			std::vector<node> nodesThatWantCommunities;
			if (totalMissingMembers < totalMissingMemberships) {
				timer.start();
				nodesThatWantCommunities.reserve(totalMissingMemberships);

				// FIXME: this depends on hash map iteration order!!!
				for (auto it : nodesWithMissingCommunities) {
					count wanted = it.second.toAssign;
					if (wanted == it.second.desired) {
						--wanted;
					}

					for (index i = 0; i < wanted; ++i) {
						nodesThatWantCommunities.push_back(it.first);
					}
				}

				const count membershipsToRemove = totalMissingMemberships - totalMissingMembers;
				while (nodesThatWantCommunities.size() > membershipsToRemove) {
					index i = drawIndex(nodesThatWantCommunities.size());
					nodesThatWantCommunities[i] = nodesThatWantCommunities.back();
					nodesThatWantCommunities.pop_back();
				}

				for (node u : nodesThatWantCommunities) {
					--nodesWithMissingCommunities[u].toAssign;
				}
				timer.stop();
				INFO("Needed ", timer.elapsedMicroseconds() * 1.0 / 1000, "ms to exclude ", nodesThatWantCommunities.size(), " memberships from assignment.");
			}

			bool overAssigned = false;
			tlx::unused(overAssigned); // needed for asserts

			// Step 3: over-assignment: sample from all nodes/memberships to get additional nodes until the numbers match
			if (totalMissingMembers > totalMissingMemberships) {
				overAssigned = true;
				timer.start();
				for (size_t i = totalMissingMemberships; i < totalMissingMembers; ++i) {
					size_t j = drawIndex(nodeMemberships.size());
					node u = nodeMemberships[j];
					// lazy deletion
					while (!nodesAlive.contains(u)) {
						nodeMemberships[j] = nodeMemberships.back();
						nodeMemberships.pop_back();
						j = drawIndex(nodeMemberships.size());
						u = nodeMemberships[j];
					}

					auto itSuccess = nodesWithMissingCommunities.insert({u, NodeMembership{desiredMemberships[u], 1}});
					if (!itSuccess.second) {
						++itSuccess.first->second.toAssign;
					}
				}
				timer.stop();
				INFO("Needed ", timer.elapsedMicroseconds() * 1.0 / 1000, "ms to sample ", totalMissingMembers - totalMissingMemberships, " additional assignments.");
			}

			timer.start();

			/**
			 * Simple bucket PQ that stores communities by desired members
			 */
			struct CommunityPQ {
				using value_type = std::pair<CommunityPtr, count>;
				// Communities sorted by desired members in decreasing order
				std::vector<value_type> communitiesByDesiredMembers;
				// Start of each bucket, i.e., range of communities with that many members
				std::vector<count> bucketStart;

				CommunityPQ(const Aux::SamplingSet<CommunityPtr> &coms) {
					if (coms.empty()) return;

					for (const CommunityPtr &com : coms) {
						const count desired = com->getDesiredNumberOfNodes() - com->getNumberOfNodes();

						if (desired >= bucketStart.size()) {
							bucketStart.resize(desired + 1);
						}

						++bucketStart[desired];
					}

					count sum = 0;
					// Reverse buckets by iterating in reverse order
					for (auto it = bucketStart.rbegin(); it != bucketStart.rend(); ++it) {
						const count tmp = *it;
						*it = sum;
						sum += tmp;
					}

					communitiesByDesiredMembers.resize(sum);

					for (const CommunityPtr &com : coms) {
						const count desired = com->getDesiredNumberOfNodes() - com->getNumberOfNodes();
						communitiesByDesiredMembers[bucketStart[desired]] = {com, desired};
						++bucketStart[desired];
					}

					for (size_t i = 0; i + 1 < bucketStart.size(); ++i) {
						bucketStart[i] = bucketStart[i + 1];
					}

					bucketStart.back() = 0;

					verify_invariants();
				}

				size_t size() const {
					return communitiesByDesiredMembers.size();
				}

				void verify_invariants() const {
#ifndef NDEBUG
					count last_desired = std::numeric_limits<count>::max();
					for (size_t i = 0; i < size(); ++i) {
						const count desired = operator[](i).second;
						assert(desired > 0);
						// Ensure that elements are in decreasing order
						if (desired != last_desired) {
							assert(desired < last_desired);
							last_desired = desired;
						}
						// assert that the item is in the right bucket
						assert(bucketStart[desired] <= i);
						assert(bucketStart[desired-1] > i);
					}
#endif
				}

                                const value_type& operator[](size_t i) const {
					return communitiesByDesiredMembers[i];
				}

				/**
				 * Decrement the counter of item at position @a pos
				 *
				 * The item is removed if the counter reaches 0.
				 * This may influence the order of elements at positions >= @a pos.
				 *
				 * @param pos The position of the item
				 */
				void decrementAt(size_t pos) {
					const count desired = communitiesByDesiredMembers[pos].second;
					assert(desired > 0);
					if (desired == 1) {
						communitiesByDesiredMembers[pos] = communitiesByDesiredMembers.back();
						communitiesByDesiredMembers.pop_back();
					} else {
						--bucketStart[desired - 1];
						index newPos = bucketStart[desired - 1];
						assert(newPos >= pos);
						using std::swap;
						swap(communitiesByDesiredMembers[pos], communitiesByDesiredMembers[newPos]);
						--communitiesByDesiredMembers[newPos].second;
					}

					verify_invariants();
				}
			};

			CommunityPQ communitiesByDesiredMembers(communitiesWantingMembers);


			std::vector<node> nodesByDesiredMemberships(nodesWithMissingCommunities.size(), none);

			{
				std::vector<count> nodesPerDesired;

				for (auto it : nodesWithMissingCommunities) {
					const count desired = it.second.desired;
					if (nodesPerDesired.size() <= desired) {
						nodesPerDesired.resize(desired + 1);
					}

					++nodesPerDesired[desired];
				}
				count sum = 0;
				// Reverse prefix sum so the actual order is reversed
				for (auto it = nodesPerDesired.rbegin(); it != nodesPerDesired.rend(); ++it) {
					const count temp = *it;
					*it = sum;
					sum += temp;
				}

				for (auto it : nodesWithMissingCommunities) {
					const count desired = it.second.desired;
					nodesByDesiredMemberships[nodesPerDesired[desired]] = it.first;
					++nodesPerDesired[desired];
				}
			}
			timer.stop();
			INFO("Needed ", timer.elapsedMicroseconds() * 1.0 / 1000, "ms to sort nodes and communities");

			// Step 4: greedy assignment
			timer.start();
			//std::vector<count> additionalMembersWanted(nodesByDesiredMemberships.size(), 0);
			count stillMissingMembers = totalMissingMembers;

			// first step: assign only nodes that actually want members
			Aux::SamplingSet<std::pair<node, Community*>, NodeCommunityHash> freshAssignments;
			// Reserve enough memory - we temporarily add additional items during switching
			freshAssignments.reserve(totalMissingMembers + 2);

			auto greedilyAssignNode =
				[&](node u, count numMembers) {
					count communitiesToFind = numMembers;

					for (index i = 0; i < communitiesByDesiredMembers.size() && communitiesToFind > 0;) {
						// Iterate from back to front
						const CommunityPtr &com = communitiesByDesiredMembers[i].first;
						if (!com->hasNode(u) && freshAssignments.insert({u, com.get()}) == 1) {
							--stillMissingMembers;
							--communitiesToFind;
							communitiesByDesiredMembers.decrementAt(i);
						} else {
							++i;
						}
					}
				};


			for (node u : nodesByDesiredMemberships) {
				greedilyAssignNode(u, nodesWithMissingCommunities[u].toAssign);
			}

			timer.stop();
			INFO("Needed ", timer.elapsedMicroseconds() * 1.0 / 1000, "ms for first greedy assignment of ", nodesByDesiredMemberships.size(), " nodes to ", communitiesWantingMembers.size(), " communities, still missing ", stillMissingMembers, " members in ", communitiesByDesiredMembers.size(), " communities.");


			// second step: if communities still want nodes, find out how many memberships are missing and add additional nodes to nodesParticipating.
			if (stillMissingMembers > 0) { // avoid log message and timer when not needed
				timer.start();
				count numRounds = 0;
				while (stillMissingMembers > 0) {
					++numRounds;
					// Step 5: if greedy assignment did not succeed, get and try additional memberships one by one:
					// if we had under-assignment, from the list of nodes whose memberships were not satisfied
					if (!nodesThatWantCommunities.empty()) {
						index j = drawIndex(nodesThatWantCommunities.size());
						node u = nodesThatWantCommunities[j];
						nodesThatWantCommunities[j] = nodesThatWantCommunities.back();
						nodesThatWantCommunities.pop_back();

						greedilyAssignNode(u, 1);
					} else {
					// if we had over-assignment, from the list of all nodes
						overAssigned = true;
						index j = drawIndex(nodeMemberships.size());
						node u = nodeMemberships[j];

						greedilyAssignNode(u, 1);
					}
				}
				timer.stop();
				INFO("Needed ", timer.elapsedMicroseconds() * 1.0 / 1000, "ms for over-assignment greedy assignment in ", numRounds, " rounds.");
			}


			// third step: randomize community assignments
			timer.start();

			const size_t numFreshAssignments = freshAssignments.size();
			assert(numFreshAssignments == totalMissingMembers);

			for (count round = 0; round < 10 * totalMissingMembers; ++round) {
				assert(numFreshAssignments == freshAssignments.size());
				std::array<node, 2> uv;
				std::array<Community*, 2> com;

				// draw node/community pairsfrom freshAssignments
				std::tie(uv[0], com[0]) = freshAssignments.at(drawIndex(numFreshAssignments));
				std::tie(uv[1], com[1]) = freshAssignments.at(drawIndex(numFreshAssignments));

				// Check if we got the same node or community twice
				if (uv[0] == uv[1] || com[0] == com[1]) continue;

				// swap both communities
				if (com[0]->hasNode(uv[1]) || com[1]->hasNode(uv[0])) continue;
				if (freshAssignments.contains({uv[0], com[1]})) continue;
				if (freshAssignments.insert({uv[1], com[0]})) {
					freshAssignments.erase({uv[0], com[0]});
					freshAssignments.erase({uv[1], com[1]});
					freshAssignments.insert({uv[0], com[1]});
				}
			}
			timer.stop();
			INFO("Needed ", timer.elapsedMicroseconds() * 1.0 / 1000, "ms for shuffling ", freshAssignments.size(), " assignments.");

			timer.start();
			// fourth step: Actually assign nodes to communities, thereby eliminating any remaining duplicates.
			for (auto it : freshAssignments) {
				const node u = it.first;
				it.second->addNode(u);
				assert(overAssigned || nodeCommunities[u].size() <= desiredMemberships[u]);
			}
			timer.stop();
			INFO("Needed ", timer.elapsedMicroseconds() * 1.0 / 1000, "ms to assign ", freshAssignments.size(), " nodes to communities");

			#ifndef NDEBUG
			for (CommunityPtr com : communities) {
				const count desired = com->getDesiredNumberOfNodes();
				assert(desired >= communitySizeSampler->getMinSize());
				const count actual = com->getNumberOfNodes();
				assert(actual == desired);
			}
			#endif
		}
	}
}
