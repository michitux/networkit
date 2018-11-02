#include "CommunityMergeEvent.h"
#include "CKBDynamicImpl.h"
#include "../../auxiliary/Random.h"
#include <tlx/unused.hpp>

namespace NetworKit {
	namespace CKBDynamicImpl {
		CommunityMergeEvent::CommunityMergeEvent(CommunityPtr communityA, CommunityPtr communityB, count targetSize, double targetEdgeProbability, count numSteps, CKBDynamicImpl& generator) : CommunityChangeEvent(generator, numSteps), communities({communityA, communityB}), targetSize(targetSize), targetEdgeProbability(targetEdgeProbability) {
			targetEdgeProbabilityPerCommunity = 1 - std::sqrt(1 - targetEdgeProbability);

			for (count c = 0; c < 2; ++c) {
				communities[c]->registerEventListener(this);

				// Get a list of nodes to be added to community c
				for (node u : communities[1-c]->getNodes()) {
					if (!communities[c]->hasNode(u)) {
						nodesToAddTo[c].insert(u);
					}
				}

				communities[c]->setAvailable(false);
			}

		}

		void CommunityMergeEvent::nextStep() {

			const count estimatedOverlap = communities[0]->getNumberOfNodes() - nodesToAddTo[1].size();

			assert(estimatedOverlap == communities[1]->getNumberOfNodes() - nodesToAddTo[0].size());
			const count estimatedMergedSize = communities[0]->getNumberOfNodes() + nodesToAddTo[0].size();
			assert(estimatedMergedSize == communities[1]->getNumberOfNodes() + nodesToAddTo[1].size());

			const count totalNodesToAdd = estimatedOverlap > targetSize ? 0 : (targetSize - estimatedOverlap) / (numSteps - currentStep);
			const count totalNodesToRemove = estimatedMergedSize < targetSize ? 0 : (estimatedMergedSize - targetSize) / (numSteps - currentStep);


			auto nonOverlappingNodes = [&]() -> count {
							   return nodesToAddTo[0].size() + nodesToAddTo[1].size();
						   };

			count nodesRemoved = 0;

			if (nonOverlappingNodes() > 0 && totalNodesToRemove > 0) {
				std::array<count, 2> nodesToRemove {nodesToAddTo[0].size(), nodesToAddTo[1].size()};
				if (nodesToRemove[0] + nodesToRemove[1] > totalNodesToRemove) {
					nodesToRemove[0] = nodesToAddTo[0].size() * totalNodesToRemove / nonOverlappingNodes();
					nodesToRemove[1] = totalNodesToRemove - nodesToRemove[0];

					if (nodesToRemove[1] > nodesToAddTo[1].size() || (nodesToAddTo[0].size() > nodesToRemove[0] && nodesToRemove[1] > 0 && Aux::Random::real() < (nodesToAddTo[0].size() *  totalNodesToRemove * 1.0 / nonOverlappingNodes()) - nodesToRemove[0])) {
						++nodesToRemove[0];
						--nodesToRemove[1];
					}
				}

				assert(nodesToRemove[0] + nodesToRemove[1] == totalNodesToRemove);
				assert(nodesToRemove[0] <= nodesToAddTo[0].size());
				assert(nodesToRemove[1] <= nodesToAddTo[1].size());

				for (count c = 0; c < 2; ++c) {
					for (count i = 0; i < nodesToRemove[c]; ++i) {
						const node u = nodesToAddTo[c].at(Aux::Random::index(nodesToAddTo[c].size()));
						nodesToAddTo[c].erase(u);

						communities[1-c]->removeNode(u);
						++nodesRemoved;
					}
				}
			}

			if (nodesRemoved < totalNodesToRemove) {
				assert(nonOverlappingNodes() == 0);
				for (; nodesRemoved < totalNodesToRemove; ++nodesRemoved) {
					assert(communities[0]->getNumberOfNodes() > 0);
					assert(communities[0]->getNumberOfNodes() == communities[1]->getNumberOfNodes());
					const node u = communities[0]->removeRandomNode();
					communities[1]->removeNode(u);
				}
			}


			for (count c = 0; c < 2; ++c) {
				// first adapt the probability so new
				// nodes get directly the right amount
				// of neighbors
				adaptProbability(communities[c], targetEdgeProbabilityPerCommunity);
			}

			count nodesAdded = 0;

			// At minimum, we need to add enough nodes such that both communities reach the minimum size.
			const count minSize = generator.communitySizeSampler->getMinSize();
			std::array<count, 2> minNodesToAdd {0,0};
			for (count c = 0; c < 2; ++c) {
				if (communities[c]->getNumberOfNodes() < minSize) {
					minNodesToAdd[c] = minSize - communities[c]->getNumberOfNodes();
				}
			}

			assert(communities[0]->getNumberOfNodes() + nodesToAddTo[0].size() == communities[1]->getNumberOfNodes() + nodesToAddTo[1].size());

			if (nonOverlappingNodes() > 0 && (totalNodesToAdd > 0 || minNodesToAdd[0] + minNodesToAdd[1] > 0)) {
				std::array<count, 2> numNodesToAdd {nodesToAddTo[0].size(), nodesToAddTo[1].size()};

				// If we have more non-overlapping nodes than nodes to be added, select some.
				// However, select only adding all non-overlapping nodes achieves the minimum size.
				if (numNodesToAdd[0] + numNodesToAdd[1] > totalNodesToAdd && communities[0]->getNumberOfNodes() + nodesToAddTo[0].size() > minSize) {
					numNodesToAdd[0] = numNodesToAdd[0] * totalNodesToAdd / nonOverlappingNodes();
					numNodesToAdd[1] = totalNodesToAdd - numNodesToAdd[0];

					// Always round the other way if the previous allocation is infeasible.
					// Apply probabilitic rounding if there is enough possibility for change.
					if (numNodesToAdd[1] > nodesToAddTo[1].size() || (nodesToAddTo[0].size() > numNodesToAdd[0] && numNodesToAdd[1] > 0 && Aux::Random::real() < nodesToAddTo[0].size() *  totalNodesToAdd *  1.0 / nonOverlappingNodes() - numNodesToAdd[0])) {
						++numNodesToAdd[0];
						--numNodesToAdd[1];
					}

					// Ensure that if possible, we are adding at least the minimum amount of nodes required to reach the minimum community size.
					for (count c = 0; c < 2; ++c) {
						if (numNodesToAdd[c] < minNodesToAdd[c]) {
							numNodesToAdd[c] = minNodesToAdd[c];
							numNodesToAdd[1-c] = std::max(totalNodesToAdd > numNodesToAdd[c] ? totalNodesToAdd - numNodesToAdd[c] : 0, minNodesToAdd[1 - c]);
						}
					}
					assert(numNodesToAdd[0] + numNodesToAdd[1] >= totalNodesToAdd);
					assert(numNodesToAdd[0] <= nodesToAddTo[0].size());
					assert(numNodesToAdd[1] <= nodesToAddTo[1].size());
				}

				for (count c = 0; c < 2; ++c) {
					for (count i = 0; i < numNodesToAdd[c]; ++i) {
						const node u = nodesToAddTo[c].at(Aux::Random::index(nodesToAddTo[c].size()));
						nodesToAddTo[c].erase(u);

						assert(communities[1-c]->hasNode(u));
						communities[c]->addNode(u);
						++nodesAdded;
					}
				}
			}

			if (nodesAdded < totalNodesToAdd || communities[0]->getNumberOfNodes() < minSize) {
				// Assert the overlap is empty now
				assert(nonOverlappingNodes() == 0);
				assert(communities[0]->getNumberOfNodes() == communities[1]->getNumberOfNodes());
				count numNodesToAdd = totalNodesToAdd - nodesAdded;
				if (communities[0]->getNumberOfNodes() + numNodesToAdd < minSize) {
					numNodesToAdd = minSize - communities[0]->getNumberOfNodes();
				}

				std::vector<node> nodes = generator.communityNodeSampler.birthCommunityNodes(numNodesToAdd, communities[0]->getNodes());
				if (nodes.size() < numNodesToAdd) {
					// Delay event completion if not enough nodes could be sampled.
					if (currentStep + 1 == numSteps) {
						++numSteps;
					}
				}

				// If we are in the really unlucky situation that we could not even get enough nodes to achieve the minimum size,
				// let both communities die and let the whole event die.
				if (communities[0]->getNumberOfNodes() + nodes.size() < minSize) {
					for (count c = 0; c < 2; ++c) {
						communities[c]->unregisterEventListener(this);
						while (communities[c]->getNumberOfNodes() > 0) {
							communities[c]->removeRandomNode();
						}
						generator.removeCommunity(communities[c]);
					}
					active = false;
					return;
				} else {
					for (node u : nodes) {
						communities[0]->addNode(u);
						communities[1]->addNode(u);
					}
				}
			}

			++currentStep;
			if (currentStep == numSteps) {
				for (count c = 0; c < 2; ++c) {
					communities[c]->unregisterEventListener(this);
				}

				active = false;
				const count oldNodes = communities[0]->getNumberOfNodes();
				tlx::unused(oldNodes);
				communities[0]->combineWith(*communities[1]);
				generator.removeCommunity(communities[1]);
				assert(communities[0]->getNumberOfNodes() == oldNodes);
				// This shouldn't change much but otherwise the community will loose half of
				// the edges in the next perturbation.
				communities[0]->changeEdgeProbability(targetEdgeProbability);
				communities[0]->setAvailable(true);
				assert(communities[0]->getNumberOfNodes() == targetSize);
			}
		}

		void CommunityMergeEvent::notifyNodeRemovedFromCommunity(node u, CommunityPtr com) {
			assert(com == communities[0] || com == communities[1]);

			for (count c = 0; c < 2; ++c) {
				nodesToAddTo[c].erase(u);
			}
		}
	}
}
