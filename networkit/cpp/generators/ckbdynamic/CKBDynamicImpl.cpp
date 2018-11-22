#include "CKBDynamicImpl.h"
#include "CommunityDeathEvent.h"
#include "CommunityBirthEvent.h"
#include "CommunitySplitEvent.h"
#include "CommunityMergeEvent.h"
#include "PowerlawCommunitySizeDistribution.h"
#include "../../auxiliary/SignalHandling.h"
#include <tlx/unused.hpp>

namespace NetworKit {
	namespace CKBDynamicImpl {
		namespace {
			/**
			 * Returns number of steps you need to wait until the next success (edge) occurs.
			 */
			count get_next_edge_distance(const double log_cp) {
				return static_cast<count>(1 + floor(log(1.0 - Aux::Random::probability()) / log_cp));
			}

		}

		void CKBDynamicImpl::addEdge(node u, node v) {
			auto e = Community::canonicalEdge(u, v);
			auto it = edgesAlive.find(e);

			if (it == edgesAlive.end()) {
				edgesAlive.insert(std::make_pair(e, 1ul));
				auto evIt = currentEdgeEvents.find(e);
				if (evIt != currentEdgeEvents.end()) {
					assert(evIt->second.type == GraphEvent::EDGE_REMOVAL);
					currentEdgeEvents.erase(evIt);
				} else {
					currentEdgeEvents.emplace(e, GraphEvent(GraphEvent::EDGE_ADDITION, e.first, e.second));
				}
			} else {
				it->second += 1;
			}
		}

		void CKBDynamicImpl::removeEdge(node u, node v) {
			auto e = Community::canonicalEdge(u, v);
			auto it = edgesAlive.find(e);
			if (it == edgesAlive.end()) {
				throw std::runtime_error("Error, removing edge that does not exist");
			}

			if (it->second > 1) {
				--it->second;
			} else {
				edgesAlive.erase(it);
				auto evIt = currentEdgeEvents.find(e);
				if (evIt != currentEdgeEvents.end()) {
					assert(evIt->second.type == GraphEvent::EDGE_ADDITION);
					currentEdgeEvents.erase(evIt);
				} else {
					currentEdgeEvents.emplace(e, GraphEvent(GraphEvent::EDGE_REMOVAL, e.first, e.second));
				}
			}
		}

		void CKBDynamicImpl::addNodeToCommunity(node u, CommunityPtr com) {
			if (com != globalCommunity) {
				nodeCommunities[u].insert(com);
				auto comNode = std::make_pair(com->getId(), u);
				auto it = currentCommunityEvents.find(comNode);
				if (it != currentCommunityEvents.end()) {
					assert(it->second.type == CommunityEvent::NODE_LEAVES_COMMUNITY);
					currentCommunityEvents.erase(it);
				} else {
					currentCommunityEvents.emplace(comNode, CommunityEvent(CommunityEvent::NODE_JOINS_COMMUNITY, u, com->getId()));
				}

				++currentCommunityMemberships;
			}
		}

		void CKBDynamicImpl::removeNodeFromCommunity(node u, CommunityPtr com) {
			if (com != globalCommunity) {
				nodeCommunities[u].erase(com);
				auto comNode = std::make_pair(com->getId(), u);
				auto it = currentCommunityEvents.find(comNode);
				if (it != currentCommunityEvents.end()) {
					assert(it->second.type == CommunityEvent::NODE_JOINS_COMMUNITY);
					currentCommunityEvents.erase(it);
				} else {
					currentCommunityEvents.emplace(comNode, CommunityEvent(CommunityEvent::NODE_LEAVES_COMMUNITY, u, com->getId()));
				}
				--currentCommunityMemberships;
			}
		}

		void CKBDynamicImpl::addCommunity(CommunityPtr com) {
			if (com->isAvailable()) {
				availableCommunities.insert(com);
			} else {
				availableCommunities.erase(com);
			}
			communities.insert(com);
		}

		void CKBDynamicImpl::removeCommunity(CommunityPtr com) {
			availableCommunities.erase(com);
			communities.erase(com);
		}

		index CKBDynamicImpl::nextCommunityId() {
			index result = maxCommunityId;
			++maxCommunityId;
			return result;
		}

		CKBDynamicImpl::CKBDynamicImpl(const CKBDynamic::param_type &params) :
			communitySizeSampler(new PowerlawCommunitySizeDistribution(params.minCommunitySize, params.maxCommunitySize, params.communitySizeExponent, params.intraCommunityEdgeProbability, params.intraCommunityEdgeExponent, params.minSplitRatio)),
			membershipDistribution(params.minCommunityMembership, params.maxCommunityMembership, params.communityMembershipExponent),
			maxCommunityId(0),
			sumOfDesiredMemberships(0),
			n(params.n),
			communityEventProbability(params.communityEventProbability),
			nodeEventProbability(params.nodeEventProbability),
			perturbationProbability(params.perturbationProbability),
			epsilon(params.epsilon),
			numTimesteps(params.numTimesteps),
			currentCommunityMemberships(0) {
			membershipDistribution.run();
		}

		std::vector<GraphEvent> CKBDynamicImpl::getGraphEvents() const {
			this->assureFinished();
			return std::move(graphEvents);
		}

		std::vector<CommunityEvent> CKBDynamicImpl::getCommunityEvents() const {
			this->assureFinished();
			return std::move(communityEvents);
		}

		void CKBDynamicImpl::generateNode() {
			node u = desiredMemberships.size();
			desiredMemberships.push_back(membershipDistribution.getDegree());
			sumOfDesiredMemberships += desiredMemberships.back();
			nodesAlive.insert(u);
			nodeCommunities.emplace_back();
			globalCommunity->addNode(u);
			graphEvents.emplace_back(GraphEvent::NODE_ADDITION, u);
		}

		void CKBDynamicImpl::eraseNode() {
			node u = nodesAlive.at(Aux::Random::index(nodesAlive.size()));
			sumOfDesiredMemberships -= desiredMemberships[u];
			desiredMemberships[u] = 0;

			while (nodeCommunities[u].size() > 0) {
				CommunityPtr com = *nodeCommunities[u].begin();
				com->removeNode(u);
			}

			assert(nodesAlive.contains(u));
			nodesAlive.erase(u);
			globalCommunity->removeNode(u);
			currentErasedNodes.push_back(u);
		}

		void CKBDynamicImpl::finishTimeStep() {
			for (const auto &it : currentEdgeEvents) {
				graphEvents.push_back(it.second);
			}
			currentEdgeEvents.clear();

			for (const auto &it : currentCommunityEvents) {
				communityEvents.push_back(it.second);
			}
			currentCommunityEvents.clear();

			// Ensure that node removals are after edge events.
			for (node u : currentErasedNodes) {
				graphEvents.emplace_back(GraphEvent::NODE_REMOVAL, u);
			}
			currentErasedNodes.clear();

			communityEvents.emplace_back(CommunityEvent::TIME_STEP);
			graphEvents.emplace_back(GraphEvent::TIME_STEP);
		}

		count CKBDynamicImpl::sampleNumSteps() const {
			return Aux::Random::integer(5, 15);
		}


		void CKBDynamicImpl::run() {
			if (hasRun) throw std::runtime_error("Error, run has already been called");

			Aux::SignalHandler handler;

			// initialization
			globalCommunity = CommunityPtr(new Community(epsilon, *this));
			communities.erase(globalCommunity);
			assert(!globalCommunity->isAvailable());

			for (node u = 0; u < n; ++u) {
				generateNode();
			}

			const count initialNumberOfNodes = nodesAlive.size();

			count sumOfDesiredMembers = 0;

			while (sumOfDesiredMembers < sumOfDesiredMemberships) {
				handler.assureRunning();
				count communitySize;
				double edgeProbability;
				std::tie(communitySize, edgeProbability) = communitySizeSampler->drawCommunity();

				CommunityPtr com(new Community(edgeProbability, *this));
				com->setDesiredNumberOfNodes(communitySize);
				com->setAvailable(true);
				sumOfDesiredMembers += communitySize;
			}

			assignNodesToCommunities();

			// Finish initial graph generation.
			finishTimeStep();

			std::binomial_distribution<count> numEventDistribution;
			double deathProbability = 0.25, birthProbability = 0.25, splitProbability = 0.25, mergeProbability = 0.25;
			tlx::unused(mergeProbability);

			for (count timestep = 0; timestep < numTimesteps; ++timestep) {
				handler.assureRunning();
				numEventDistribution.param(std::binomial_distribution<count>::param_type(communities.size(), communityEventProbability));
				const count numCommunityEvents = numEventDistribution(Aux::Random::getURNG());

				numEventDistribution.param(std::binomial_distribution<count>::param_type(nodesAlive.size(), nodeEventProbability));
				const count numNodeEvents = numEventDistribution(Aux::Random::getURNG());

				INFO("Timestep ", timestep, " generating ", numCommunityEvents, " community events and ", numNodeEvents, " node events");

				for (count i = 0; i < numCommunityEvents; ++i) {
					handler.assureRunning();
					count numSteps = sampleNumSteps();
					double r = Aux::Random::real();
					if (r < birthProbability) {
						// generate new community
						currentEvents.emplace_back(new CommunityBirthEvent(numSteps, *this));
					} else if (r < birthProbability + deathProbability) {
						// let a community die
						if (availableCommunities.size() > 0) {
							CommunityPtr com = availableCommunities.at(Aux::Random::index(availableCommunities.size()));
							count coreSize = std::max<count>(0.1 * com->getNumberOfNodes(), communitySizeSampler->getMinSize());
							currentEvents.emplace_back(new CommunityDeathEvent(com, coreSize, numSteps, *this));
							assert(!com->isAvailable());
						} else {
							WARN("No community available for death event.");
						}
					} else if (r < birthProbability + deathProbability + splitProbability) {
						// Split a community
						if (availableCommunities.size() > 0) {
							CommunityPtr com = availableCommunities.at(Aux::Random::index(availableCommunities.size()));
							auto comSizeProbA = communitySizeSampler->drawCommunity();
							auto comSizeProbB = communitySizeSampler->drawCommunity();
							currentEvents.emplace_back(new CommunitySplitEvent(com, comSizeProbA.first, comSizeProbA.second, comSizeProbB.first, comSizeProbB.second, numSteps, *this));
							assert(!com->isAvailable());
						} else {
							WARN("No community available for splitting.");
						}
					} else {
						// merge two communities
						if (availableCommunities.size() > 1) {
							index ia = Aux::Random::integer(availableCommunities.size() - 1);
							index ib = Aux::Random::integer(1, availableCommunities.size() - 1);
							if (ia == ib) {
								ib = 0;
							}

							CommunityPtr comA = availableCommunities.at(ia);
							CommunityPtr comB = availableCommunities.at(ib);

							count targetSize;
							double targetEdgeProbability;
							std::tie(targetSize, targetEdgeProbability) = communitySizeSampler->drawCommunity();
							currentEvents.emplace_back(new CommunityMergeEvent(comA, comB, targetSize, targetEdgeProbability, numSteps, *this));
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
				numEventDistribution.param(std::binomial_distribution<count>::param_type(numNodeEvents, nodeBirthProbability));
				const count nodesBorn = numEventDistribution(Aux::Random::getURNG());

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
					const double x = sumOfDesiredMemberships * 1.0 / currentCommunityMemberships;
					birthProbability = 0.5 * x / (1 + x);
					deathProbability = 0.5 - birthProbability;
					INFO("Current memberships: ", currentCommunityMemberships, " desired: ", sumOfDesiredMemberships, " number of communities: ", communities.size(), " available: ", availableCommunities.size(), " active events ", currentEvents.size());
					INFO("Current nodes ", nodesAlive.size(), " current edges: ", edgesAlive.size(), " total graph events ", graphEvents.size(), " total community events ", communityEvents.size());
				}

				finishTimeStep();

			}

			hasRun = true;
		}

		void CKBDynamicImpl::assignNodesToCommunities() {
			count totalMissingMembers = 0;

			std::vector<CommunityPtr> communitiesWithMissingMembers;

			for (CommunityPtr com : communities) {
				const count desired = com->getDesiredNumberOfNodes();
				const count actual = com->getNumberOfNodes();
				assert(actual <= desired);

				if (actual < desired) {
					communitiesWithMissingMembers.push_back(com);
					totalMissingMembers += desired - actual;
				}
			}

			if (totalMissingMembers == 0) return;

			count totalMissingMemberships = 0;

			for (node u : nodesAlive) {
				const count desired = desiredMemberships[u];
				const count actual = nodeCommunities[u].size();

				if (desired > actual) {
					totalMissingMemberships += desired - actual;
				}
			}

			// If totalMissingMemberships > totalMissingMembers, a) find nodes that have too many memberships and remove those nodes from some of them and b) create empty communities of size 1 that are removed afterwards.
			// If totalMissingMembers > totalMissingMemberships, find nodes to which we can assign additional memberships, i.e., nodes that are not much over their allocated memberships.

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
				std::vector<std::pair<node, CommunityPtr>> nodesWithAdditionalMemberships;
				std::vector<CommunityPtr> ac;
				for (node u : nodesAlive) {
					if (desiredMemberships[u] < nodeCommunities[u].size()) {
						count numTooMany = nodeCommunities[u].size() - desiredMemberships[u];
						for (const CommunityPtr &com : nodeCommunities[u]) {
							if (com->isAvailable()) {
								ac.push_back(com);
							}
						}

						while (ac.size() > numTooMany) {
							const index i = Aux::Random::index(ac.size());
							std::swap(ac[i], ac.back());
							ac.pop_back();
						}

						for (CommunityPtr com : ac) {
							nodesWithAdditionalMemberships.emplace_back(u, com);
						}

						ac.clear();
					}
				}

				std::shuffle(nodesWithAdditionalMemberships.begin(), nodesWithAdditionalMemberships.end(), Aux::Random::getURNG());

				while (!nodesWithAdditionalMemberships.empty() && totalMissingMembers < totalMissingMemberships) {
					const std::pair<node, CommunityPtr> u_com = nodesWithAdditionalMemberships.back();
					nodesWithAdditionalMemberships.pop_back();

					bool alreadyListed = (u_com.second->getDesiredNumberOfNodes() > u_com.second->getNumberOfNodes());

					u_com.second->removeNode(u_com.first);
					++totalMissingMembers;

					if (!alreadyListed) {
						communitiesWithMissingMembers.push_back(u_com.second);
					}
				}
			}


			std::vector<std::pair<CommunityPtr, count>> communitiesByDesiredMembers;

			{
				std::vector<count> numCommunitiesWithDesired;

				for (CommunityPtr com : communitiesWithMissingMembers) {
					const count desired = com->getDesiredNumberOfNodes();

					if (desired >= numCommunitiesWithDesired.size()) {
						numCommunitiesWithDesired.resize(desired + 1);
					}

					++numCommunitiesWithDesired[desired];
				}

				count sum = 0;
				for (auto it = numCommunitiesWithDesired.begin(); it != numCommunitiesWithDesired.end(); ++it) {
					const count tmp = *it;
					*it = sum;
					sum += tmp;
				}

				communitiesByDesiredMembers.resize(sum);

				for (CommunityPtr com : communitiesWithMissingMembers) {
					const count desired = com->getDesiredNumberOfNodes();
					const count actual = com->getNumberOfNodes();
					communitiesByDesiredMembers[numCommunitiesWithDesired[desired]] = {com, desired - actual};
					++numCommunitiesWithDesired[desired];
				}
			}


			std::vector<node> nodesByDesiredMemberships(nodesAlive.size(), none);

			{
				std::vector<count> nodesPerDesired;

				for (node u : nodesAlive) {
					const count desired = desiredMemberships[u];
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

				for (node u : nodesAlive) {
					const count desired = desiredMemberships[u];
					nodesByDesiredMemberships[nodesPerDesired[desired]] = u;
					++nodesPerDesired[desired];
				}
			}


			std::vector<node> nodesParticipating;
			std::vector<bool> nodeIsParticipating(nodesByDesiredMemberships.size(), false);
			count stillMissingMembers = totalMissingMembers;

			// first step: assign only nodes that actually want members
			std::vector<Aux::SamplingSet<CommunityPtr>> freshAssignments(nodesByDesiredMemberships.size());
			// How many communities shall be assigned in the greedy assignment step.
			std::vector<count> numInitialCommunitiesToAssign(nodesByDesiredMemberships.size(), 0);

			for (node lu = 0; lu < nodesByDesiredMemberships.size(); ++lu) {
				const node u = nodesByDesiredMemberships[lu];
				if (desiredMemberships[u] > nodeCommunities[u].size()) {
					nodesParticipating.push_back(lu);
					nodeIsParticipating[lu] = true;
					numInitialCommunitiesToAssign[lu] = desiredMemberships[u] - nodeCommunities[u].size();
				}
			}

			double overAssignment = 0;

			auto greedilyAssignNodes =
				[&]() {
					for (node lu = 0; lu < nodesByDesiredMemberships.size(); ++lu) {
						const node u = nodesByDesiredMemberships[lu];
						if (nodeIsParticipating[lu]) {
							count communitiesToFind = numInitialCommunitiesToAssign[lu] - freshAssignments[lu].size();
							index last_empty = communitiesByDesiredMembers.size();
							for (index i = 0; i < communitiesByDesiredMembers.size() && communitiesToFind > 0; ++i) {
								const index ci = communitiesByDesiredMembers.size() - i - 1;
								const CommunityPtr &com = communitiesByDesiredMembers[ci].first;
								count &missing = communitiesByDesiredMembers[ci].second;
								if (missing > 0 && !com->hasNode(u)) {
									if (freshAssignments[lu].insert(com) == 1) {
										--missing;
										--stillMissingMembers;
										--communitiesToFind;
										if (communitiesToFind == 0) {
											break;
										}
									}
								}

								if (missing == 0) {
									last_empty = ci;
								}
							}

							if (last_empty < communitiesByDesiredMembers.size()) {
								index wi = last_empty;
								for (index ri = last_empty; ri < communitiesByDesiredMembers.size(); ++ri) {
									if (communitiesByDesiredMembers[ri].second > 0) {
										communitiesByDesiredMembers[wi] = communitiesByDesiredMembers[ri];
										++wi;
									}
								}

								while (communitiesByDesiredMembers.size() > wi) {
									communitiesByDesiredMembers.pop_back();
								}

							}

							assert(freshAssignments[lu].size() + nodeCommunities[u].size() <= desiredMemberships[u] || overAssignment > 0);
						}
					}
				};

			greedilyAssignNodes();


			// second step: if communities still want nodes, find out how many memberships are missing and add additional nodes to nodesParticipating.
			// FIXME: should we initialize overAssignment differently and possibly already in the very first assignment step?
			while (stillMissingMembers > 0) {
				overAssignment += std::max(0.01, stillMissingMembers * 1.0 / sumOfDesiredMemberships);

				for (node lu = 0; lu < nodesByDesiredMemberships.size(); ++lu) {
					const node u = nodesByDesiredMemberships[lu];

					// Desired memberships with over assignment
					double fdwo = desiredMemberships[u] * (1.0 + overAssignment);
					count dwo = desiredMemberships[u] * (1.0 + overAssignment);
					if (Aux::Random::real() < fdwo - dwo) {
						++dwo;
					}

					if (dwo > nodeCommunities[u].size() + numInitialCommunitiesToAssign[lu]) {
						if (!nodeIsParticipating[lu]) {
							nodeIsParticipating[lu] = true;
							nodesParticipating.push_back(lu);
						}

						numInitialCommunitiesToAssign[lu] = dwo - nodeCommunities[u].size();
					}
				}

				greedilyAssignNodes();
			}

			// third step: randomize community assignments of nodesParticipating, balance assignments if there were over-assignments.

			for (count round = 0; round < 100; ++round) {
				std::shuffle(nodesParticipating.begin(), nodesParticipating.end(), Aux::Random::getURNG());

				for (count i = 0; i + 1 < nodesParticipating.size(); i += 2) {
					const std::array<node, 2> ln {nodesParticipating[i], nodesParticipating[i+1]};
					const std::array<node, 2> uv {nodesByDesiredMemberships[ln[0]], nodesByDesiredMemberships[ln[1]]};

					std::vector<CommunityPtr> communitiesToShuffle;
					for (count j = 0; j < 2; ++j) {
						for (CommunityPtr com : freshAssignments[ln[j]]) {
							if (!freshAssignments[ln[1-j]].contains(com) && !com->hasNode(uv[1-j])) {
								communitiesToShuffle.push_back(com);
								freshAssignments[ln[j]].erase(com);
							}
						}
					}

					if (communitiesToShuffle.size() > 0) {
						std::shuffle(communitiesToShuffle.begin(), communitiesToShuffle.end(), Aux::Random::getURNG());

						// Assure that, if possible, no node has no community at the end
						for (count j = 0; j < 2 && !communitiesToShuffle.empty(); ++j) {
							if (freshAssignments[ln[j]].empty() && nodeCommunities[uv[j]].empty()) {
								freshAssignments[ln[j]].insert(communitiesToShuffle.back());
								communitiesToShuffle.pop_back();
							}
						}

						while (!communitiesToShuffle.empty()) {
							// Calculate the percentage of the desired memberships we would get if we assigned the community
							// to the first/second node.
							std::array<double, 2> possibleOverAssignment;

							for (index j = 0; j < 2; ++j) {
								possibleOverAssignment[j] = (freshAssignments[ln[j]].size() + nodeCommunities[uv[j]].size() + 1.0) * 1.0 / desiredMemberships[uv[j]];
							}

							index toAssign = 0;
							if (possibleOverAssignment[0] > possibleOverAssignment[1]) {
								toAssign = 1;
							}

							freshAssignments[ln[toAssign]].insert(communitiesToShuffle.back());
							communitiesToShuffle.pop_back();
						}

					}
				}
			}


			// fourth step: Actually assign nodes to communities, thereby eliminating any remaining duplicates.
			for (node lu : nodesParticipating) {
				const node u = nodesByDesiredMemberships[lu];
				for (CommunityPtr com : freshAssignments[lu]) {
					com->addNode(u);
				}
				assert(nodeCommunities[u].size() <= desiredMemberships[u] || overAssignment > 0);
			}
		}
	}
}
