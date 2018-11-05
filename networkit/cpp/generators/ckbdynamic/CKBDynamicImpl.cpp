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
				communityNodeSampler.assignCommunity(u);
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
				communityNodeSampler.leaveCommunity(u);
				--currentCommunityMemberships;
			}
		}

		void CKBDynamicImpl::addCommunity(CommunityPtr com) {
			if (com->isAvailable()) {
				// If community is too small, remove community again!!
				if (com->getNumberOfNodes() < communitySizeSampler->getMinSize()) {
					INFO("community has only ", com->getNumberOfNodes(), " nodes, destroying.");
					currentEvents.emplace_back(new CommunityDeathEvent(com, 0, 1, *this));
				} else {
					availableCommunities.insert(com);
				}
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
			communityNodeSampler(0, params.minCommunityMembership, params.maxCommunityMembership, params.communityMembershipExponent),
			communitySizeSampler(new PowerlawCommunitySizeDistribution(params.minCommunitySize, params.maxCommunitySize, params.communitySizeExponent, params.intraCommunityEdgeProbability, params.intraCommunityEdgeExponent, params.minSplitRatio)),
			maxCommunityId(0),
			n(params.n),
			communityEventProbability(params.communityEventProbability),
			nodeEventProbability(params.nodeEventProbability),
			perturbationProbability(params.perturbationProbability),
			epsilon(params.epsilon),
			numTimesteps(params.numTimesteps),
			currentCommunityMemberships(0) {
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
			node u = communityNodeSampler.addNode();
			for (node v = u; v < communityNodeSampler.getNumberOfNodes(); ++v) {
				nodesAlive.insert(v);
				nodeCommunities.emplace_back();
				globalCommunity->addNode(v);
				graphEvents.emplace_back(GraphEvent::NODE_ADDITION, v);
			}
		}

		void CKBDynamicImpl::eraseNode() {
			std::vector<node> removedNodes;
			removedNodes.push_back(nodesAlive.at(Aux::Random::index(nodesAlive.size())));
			node v = communityNodeSampler.removeNode(removedNodes.back());
			// Also remove the additionally removed node if there was any.
			if (v != none) removedNodes.push_back(v);

			for (node u : removedNodes) {
				while (nodeCommunities[u].size() > 0) {
					CommunityPtr com = *nodeCommunities[u].begin();
					com->removeNode(u);
					// if a community becomes too small, erase it
					if (com->isAvailable() && com->getNumberOfNodes() < communitySizeSampler->getMinSize()) {
						INFO("Available community has only ", com->getNumberOfNodes(), " nodes, destroying.");
						currentEvents.emplace_back(new CommunityDeathEvent(com, 0, 1, *this));
					}
				}

				assert(nodesAlive.contains(u));
				nodesAlive.erase(u);
				globalCommunity->removeNode(u);
				currentErasedNodes.push_back(u);
			}
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

			while (currentCommunityMemberships < communityNodeSampler.getSumOfDesiredMemberships()) {
				handler.assureRunning();
				count communitySize;
				double edgeProbability;
				std::tie(communitySize, edgeProbability) = communitySizeSampler->drawCommunity();
				std::vector<node> comNodes(communityNodeSampler.birthCommunityNodes(communitySize));
				CommunityPtr com(new Community(edgeProbability, *this));
				for (node u : comNodes) {
					com->addNode(u);
				}
				com->setAvailable(true);
			}

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

					const double sqrtPerturbationProbability = std::sqrt(perturbationProbability);

					const double log_perturb = std::log(1.0 - sqrtPerturbationProbability);

					for (count ci = get_next_edge_distance(log_perturb) - 1; ci < communities.size(); ci += get_next_edge_distance(log_perturb)) {
						handler.assureRunning();
						communities.at(ci)->perturbEdges(sqrtPerturbationProbability);
					}
				}

				// adjust event probabilities
				{
					const double x = communityNodeSampler.getSumOfDesiredMemberships() * 1.0 / currentCommunityMemberships;
					birthProbability = 0.5 * x / (1 + x);
					deathProbability = 0.5 - birthProbability;
					INFO("Current memberships: ", currentCommunityMemberships, " desired: ", communityNodeSampler.getSumOfDesiredMemberships(), " number of communities: ", communities.size(), " available: ", availableCommunities.size(), " active events ", currentEvents.size());
					INFO("Current nodes ", nodesAlive.size(), " current edges: ", edgesAlive.size(), " total graph events ", graphEvents.size(), " total community events ", communityEvents.size());
				}

				finishTimeStep();

			}

			hasRun = true;
		}

	}
}
