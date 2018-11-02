#include "CommunitySplitEvent.h"

#include "CKBDynamicImpl.h"

namespace NetworKit {
	namespace CKBDynamicImpl {
		CommunitySplitEvent::CommunitySplitEvent(CommunityPtr community, count targetSizeA, double targetEdgeProbabilityA, count targetSizeB, double targetEdgeProbabilityB, count numSteps, CKBDynamicImpl& generator) : CommunityChangeEvent(generator, numSteps), targetSize({targetSizeA, targetSizeB}), targetEdgeProbability({targetEdgeProbabilityA, targetEdgeProbabilityB}), communities({community, CommunityPtr(new Community(*community))}) {
			std::vector<node> nodes;
			nodes.reserve(communities[0]->getNumberOfNodes());
			for (node u : communities[0]->getNodes()) {
				nodes.push_back(u);
			}

			std::shuffle(nodes.begin(), nodes.end(), Aux::Random::getURNG());

			while (targetSizeA + targetSizeB < nodes.size()) {
				// sample a random set of nodes to remove from both
				node u = nodes.back();
				nodes.pop_back();
				nodesToRemove[0].insert(u);
				nodesToRemove[1].insert(u);
			}

			// Calculate how many nodes need to be removed from the two communities.
			// Due to the previous loop, both numbers sum up to at least nodes.size().
			const count combinedTargetSize = targetSizeA + targetSizeB;
			assert(combinedTargetSize >= nodes.size());
			double fractionA = targetSizeA * 1.0 / combinedTargetSize;
			std::array<count, 2> numNodesToRemove {targetSizeB * nodes.size() / combinedTargetSize, 0};
			numNodesToRemove[1] = nodes.size() - numNodesToRemove[0];

			if (nodes.size() - numNodesToRemove[0] > targetSizeA) {
				numNodesToRemove[0] = nodes.size() - targetSizeA;
				numNodesToRemove[1] = nodes.size() - numNodesToRemove[0];
			} else if (nodes.size() - numNodesToRemove[1] > targetSizeB) {
				numNodesToRemove[1] = nodes.size() - targetSizeB;
				numNodesToRemove[0] = nodes.size() - numNodesToRemove[1];
			} else if (nodes.size() - numNodesToRemove[1] < targetSizeB && numNodesToRemove[1] > 0 && Aux::Random::real() < fractionA * nodes.size() - numNodesToRemove[1]) {
				++numNodesToRemove[0];
				--numNodesToRemove[1];
			}

			assert(nodes.size() - numNodesToRemove[0] <= targetSizeA);
			assert(nodes.size() - numNodesToRemove[1] <= targetSizeB);
			assert(numNodesToRemove[0] + numNodesToRemove[1] == nodes.size());

			for (count com = 0; com < 2; ++com) {
				for (count i = 0; i < numNodesToRemove[com]; ++i) {
					assert(!nodes.empty());
					node u = nodes.back();
					nodes.pop_back();
					nodesToRemove[com].insert(u);
				}
			}

			assert(nodes.empty());

			assert(communities[0]->getNumberOfNodes() - nodesToRemove[0].size() <= targetSize[0]);
			assert(communities[1]->getNumberOfNodes() - nodesToRemove[1].size() <= targetSize[1]);

			communities[0]->setAvailable(false);
			communities[1]->setAvailable(false);
			communities[0]->registerEventListener(this);
			communities[1]->registerEventListener(this);
		}

		void CommunitySplitEvent::nextStep() {
			assert(communities[0]->getNumberOfNodes() - nodesToRemove[0].size() <= targetSize[0]);
			assert(communities[1]->getNumberOfNodes() - nodesToRemove[1].size() <= targetSize[1]);

			for (count com = 0; com < 2; ++com) {
				if (communities[com]->getNumberOfNodes() == 0) continue;

				assert(communities[com]->getNumberOfNodes() - nodesToRemove[com].size() <= targetSize[com]);

				const count numNodesToRemove = nodesToRemove[com].size() / (numSteps - currentStep);
				for (count i = 0; i < numNodesToRemove; ++i) {
					const node u = nodesToRemove[com].at(Aux::Random::index(nodesToRemove[com].size()));
					nodesToRemove[com].erase(u);

					assert(generator.hasNode(u));
					communities[com]->removeNode(u);
				}

				adaptProbability(communities[com], targetEdgeProbability[com]);

				// Add nodes to achieve target size, and at least the minimum size in the current step
				if (targetSize[com] > (communities[com]->getNumberOfNodes() - nodesToRemove[com].size())) {
					count numNodesToAdd = (targetSize[com] + nodesToRemove[com].size() - communities[com]->getNumberOfNodes()) / (numSteps - currentStep);
					if (communities[com]->getNumberOfNodes() + numNodesToAdd < generator.communitySizeSampler->getMinSize()) {
						numNodesToAdd = generator.communitySizeSampler->getMinSize() - communities[com]->getNumberOfNodes();
					}

					if (numNodesToAdd > 0) {
						assert(communities[com]->getNumberOfNodes() + numNodesToAdd - nodesToRemove[com].size() <= targetSize[com]);
						const std::vector<node> nodesToAdd = generator.communityNodeSampler.birthCommunityNodes(numNodesToAdd, communities[com]->getNodes());

						if (nodesToAdd.size() < numNodesToAdd) {
							// delay event if we are in the last step
							if (currentStep + 1 == numSteps) {
								++numSteps;
							}
						}

						// After this step, the community will be too small to survive.
						// Give up this community!
						if (communities[com]->getNumberOfNodes() + nodesToAdd.size() < generator.communitySizeSampler->getMinSize()) {
							while (communities[com]->getNumberOfNodes() > 0) {
								communities[com]->removeRandomNode();
							}

							generator.removeCommunity(communities[com]);
						} else {
							for (node u : nodesToAdd) {
								communities[com]->addNode(u);
							}

                            assert(communities[com]->getNumberOfNodes() - nodesToRemove[com].size() <= targetSize[com]);
						}
					}
				}
			}

			++currentStep;
			if (currentStep == numSteps) {
				active = false;
				for (count com = 0; com < 2; ++com) {
					communities[com]->unregisterEventListener(this);

					// Do not mark community as available if it has died due to getting too small.
					if (communities[com]->getNumberOfNodes() > 0) {
						assert(communities[com]->getNumberOfNodes() == targetSize[com]);
						communities[com]->setAvailable(true);
					}
				}
			}
		}

		void CommunitySplitEvent::notifyNodeRemovedFromCommunity(node u, CommunityPtr com) {
			// Remove node from the nodes to be removed of the respective community.
			if (com == communities[0]) {
				nodesToRemove[0].erase(u);
			} else {
				assert(com == communities[1]);
				nodesToRemove[1].erase(u);
			}
		}

	}
}
