#include "CommunityMergeEvent.h"
#include "CKBDynamicImpl.h"
#include "../../auxiliary/Random.h"
#include <tlx/unused.hpp>

namespace NetworKit {
	namespace CKBDynamicImpl {
		CommunityMergeEvent::CommunityMergeEvent(CommunityPtr communityA, CommunityPtr communityB, count targetSize, count numSteps, CKBDynamicImpl& generator) : CommunityChangeEvent(generator, numSteps), communities({communityA, communityB}), targetSize(targetSize), communitiesMerged(false) {
			//targetEdgeProbabilityPerCommunity = 1 - std::sqrt(1 - targetEdgeProbability);

			for (count c = 0; c < 2; ++c) {
				communities[c]->setCurrentEvent(this);

				// Get a list of nodes to be added to community c
				for (node u : communities[1-c]->getNodes()) {
					if (!communities[c]->hasNode(u)) {
						nodesToAddTo[c].insert(u);
					}
				}
			}

		}

		void CommunityMergeEvent::nextStep() {

			const count estimatedOverlap = communities[0]->getNumberOfNodes() - nodesToAddTo[1].size();

			assert(communitiesMerged || estimatedOverlap == communities[1]->getNumberOfNodes() - nodesToAddTo[0].size());
			const count estimatedMergedSize = communities[0]->getNumberOfNodes() + nodesToAddTo[0].size();
			assert(communitiesMerged || estimatedMergedSize == communities[1]->getNumberOfNodes() + nodesToAddTo[1].size());

			const count totalNodesToAdd = estimatedOverlap > targetSize ? 0 : std::ceil(static_cast<double>(targetSize - estimatedOverlap) / (numSteps - currentStep));
			const count totalNodesToRemove = estimatedMergedSize < targetSize ? 0 : std::ceil(static_cast<double>(estimatedMergedSize - targetSize) / (numSteps - currentStep));


			auto nonOverlappingNodes = [&]() -> count {
							   return nodesToAddTo[0].size() + nodesToAddTo[1].size();
						   };

			count nodesRemoved = 0;

			// First: remove nodes that are not in the overlap if nodes need to be removed
			if (nonOverlappingNodes() > 0 && totalNodesToRemove > 0) {
				std::array<count, 2> nodesToRemove {nodesToAddTo[0].size(), nodesToAddTo[1].size()};
				if (nodesToRemove[0] + nodesToRemove[1] > totalNodesToRemove) {
					// Distribute nodes to remove proportionally to the non-overlapping parts
					nodesToRemove[0] = nodesToAddTo[0].size() * totalNodesToRemove / nonOverlappingNodes();
					nodesToRemove[1] = totalNodesToRemove - nodesToRemove[0];

					if (nodesToRemove[1] > nodesToAddTo[1].size() || (nodesToAddTo[0].size() > nodesToRemove[0] && nodesToRemove[1] > 0 && Aux::Random::real() < (nodesToAddTo[0].size() *  totalNodesToRemove * 1.0 / nonOverlappingNodes()) - nodesToRemove[0])) {
						++nodesToRemove[0];
						--nodesToRemove[1];
					}

					assert(nodesToRemove[0] + nodesToRemove[1] == totalNodesToRemove);
				}

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

			// Merge communities if we have no nodes in the overlap left
			if (nonOverlappingNodes() == 0) {
				mergeCommunities();
			}

			// If we still need to remove nodes, we are now after the merge.
			if (nodesRemoved < totalNodesToRemove) {
				assert(nonOverlappingNodes() == 0);
				for (; nodesRemoved < totalNodesToRemove; ++nodesRemoved) {
					assert(communities[0]->getNumberOfNodes() > 0);
					// communities[1] doesn't exist anymore
					communities[0]->removeRandomNode();
				}
			}


			// At minimum, we need to add enough nodes such that both communities reach the minimum size.
			const count minSize = generator.communitySizeSampler->getMinSize();

			if (communitiesMerged) {
				assert(nonOverlappingNodes() == 0);

				communities[0]->setDesiredNumberOfNodes(std::max(minSize, communities[0]->getNumberOfNodes() + totalNodesToAdd));
			} else {
				assert(communities[0]->getNumberOfNodes() + nodesToAddTo[0].size() == communities[1]->getNumberOfNodes() + nodesToAddTo[1].size());

				// Check if we need to fully merge the communities to achieve the minimum size or the desired growth of the overlap
				if (communities[0]->getNumberOfNodes() + nodesToAddTo[0].size() <= minSize || totalNodesToAdd >= nodesToAddTo[0].size() + nodesToAddTo[1].size()) {
					// Add all remaining non-overlapping nodes and perform merge immediately, might need to further grow the community.
					count sizeAfterStep = communities[0]->getNumberOfNodes() + nodesToAddTo[0].size();
					if (totalNodesToAdd > nodesToAddTo[0].size() + nodesToAddTo[1].size()) {
						sizeAfterStep += (totalNodesToAdd - (nodesToAddTo[0].size() + nodesToAddTo[1].size()));
					}

					if (sizeAfterStep < minSize) {
						sizeAfterStep = minSize;
					}

					assert(sizeAfterStep <= targetSize);

					communities[0]->setDesiredNumberOfNodes(sizeAfterStep);

					while (!nodesToAddTo[0].empty()) {
						const node u = nodesToAddTo[0].at(0);

						assert(communities[1]->hasNode(u));
						communities[0]->addNode(u);
						nodesToAddTo[0].erase(u);
					}

					// Discard the second community, there is no need to add those nodes to the second community.
					// The first community will simply replace it.
					nodesToAddTo[1].clear();

					mergeCommunities();
				} else {
					// Adding nodes from nodesToAddTo will be enough, no merge will be performed and no additional nodes will be needed

					std::array<count, 2> numNodesToAdd {0, 0};
					// First step: ensure that we add enough nodes to achieve the minimum size
					for (count c = 0; c < 2; ++c) {
						if (communities[c]->getNumberOfNodes() < minSize) {
							numNodesToAdd[c] = minSize - communities[c]->getNumberOfNodes();
						}
					}

					// Second step: add more nodes if we wanted to add more
					if (numNodesToAdd[0] + numNodesToAdd[1] < totalNodesToAdd) {
						// How many nodes we need to add additionally
						count additionalNodes = totalNodesToAdd - numNodesToAdd[0] - numNodesToAdd[1];
						// How many nodes are available
						std::array<count, 2> availableNodes {
							nodesToAddTo[0].size() - numNodesToAdd[0],
							nodesToAddTo[1].size() - numNodesToAdd[1]
						};
						count allAvailableNodes = availableNodes[0] + availableNodes[1];
						assert(additionalNodes <= allAvailableNodes);

						// What fraction of available nodes we should add
						double fractionToAdd = additionalNodes * 1.0 / allAvailableNodes;

						// Get this fraction, but round down
						std::array<count, 2> additionalNodesPerCommunity {
							// Use integer arithmetic to ensure that there are no rounding errors
							availableNodes[0] * additionalNodes / allAvailableNodes,
							availableNodes[1] * additionalNodes / allAvailableNodes
						};

						// There might be one more node missing, add it
						if (additionalNodesPerCommunity[0] + additionalNodesPerCommunity[1] < additionalNodes) {
							// probabilistic rounding with the fraction to add
							if (Aux::Random::real() < availableNodes[0] * fractionToAdd - additionalNodesPerCommunity[0]) {
								++additionalNodesPerCommunity[0];
							} else {
								++additionalNodesPerCommunity[1];
							}
						}

						assert(additionalNodes == additionalNodesPerCommunity[0] + additionalNodesPerCommunity[1]);

						numNodesToAdd[0] += additionalNodesPerCommunity[0];
						numNodesToAdd[1] += additionalNodesPerCommunity[1];
					}

					// Third step: add nodes to the community
					for (count c = 0; c < 2; ++c) {
						assert(numNodesToAdd[c] <= nodesToAddTo[c].size());
						// Adjust the size and edge probability before adding nodes
						communities[c]->setDesiredNumberOfNodes(communities[c]->getNumberOfNodes() + numNodesToAdd[c]);

						// Actually add as many nodes as desired to increase the overlap
						for (count i = 0; i < numNodesToAdd[c]; ++i) {
							const node u = nodesToAddTo[c].at(Aux::Random::index(nodesToAddTo[c].size()));

							assert(communities[1-c]->hasNode(u));
							communities[c]->addNode(u);
							nodesToAddTo[c].erase(u);
						}
					}
				}
			} // if (communitiesMerge)


			++currentStep;
			if (currentStep == numSteps) {
				mergeCommunities();

				assert(communities[0]->getDesiredNumberOfNodes() == targetSize);
				//communities[0]->changeEdgeProbability(targetEdgeProbability);
				communities[0]->setCurrentEvent(nullptr);
				active = false;
			}
		}

		void CommunityMergeEvent::mergeCommunities() {
			if (communitiesMerged) return;

			communities[1]->setCurrentEvent(nullptr);
			while (communities[1]->getNumberOfNodes() > 0) {
				communities[1]->removeRandomNode();
			}

			generator.removeCommunity(communities[1]);

			communitiesMerged = true;
		}

		void CommunityMergeEvent::notifyNodeRemovedFromCommunity(node u, CommunityPtr com) {
			tlx::unused(com);
			assert(com == communities[0] || com == communities[1]);

			for (count c = 0; c < 2; ++c) {
				nodesToAddTo[c].erase(u);
			}
		}

		void CommunityMergeEvent::notifyNodeAddedToCommunity(node u, CommunityPtr com) {
			tlx::unused(u);
			tlx::unused(com);
			if (!communitiesMerged) {
				if (com == communities[0]) {
					assert(nodesToAddTo[0].contains(u));
				} else {
					assert(nodesToAddTo[1].contains(u));
				}
			}
		}

		bool CommunityMergeEvent::canRemoveNode() const {
			return communitiesMerged;
		}
	}
}
