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
			double fractionA = targetSizeA * 1.0 / combinedTargetSize;
			std::array<count, 2> numNodesToRemove {targetSizeA * nodes.size() / combinedTargetSize, 0};
			numNodesToRemove[1] = nodes.size() - numNodesToRemove[0];

			if (nodes.size() - numNodesToRemove[0] > targetSizeA) {
				numNodesToRemove[0] = nodes.size() - targetSizeA;
				numNodesToRemove[1] = nodes.size() - numNodesToRemove[0];
			} else if (nodes.size() - numNodesToRemove[1] > targetSizeB) {
				numNodesToRemove[1] = nodes.size() - targetSizeB;
				numNodesToRemove[0] = nodes.size() - numNodesToRemove[1];
			} else if (nodes.size() - numNodesToRemove[1] < targetSizeB && numNodesToRemove[1] > 0 && Aux::Random::real() < fractionA * nodes.size() - numNodesToRemove[0]) {
				++numNodesToRemove[0];
				--numNodesToRemove[1];
			}

			assert(nodes.size() - numNodesToRemove[0] <= targetSizeA);
			assert(nodes.size() - numNodesToRemove[1] <= targetSizeB);
			assert(numNodesToRemove[0] + numNodesToRemove[1] == nodes.size());

			for (count com = 0; com < 2; ++com) {
				for (count i = 0; i < numNodesToRemove[com]; ++i) {
					node u = nodes.back();
					nodes.pop_back();
					nodesToRemove[com].insert(u);
				}
			}

			communities[0]->setAvailable(false);
			communities[1]->setAvailable(false);
			communities[0]->registerEventListener(this);
			communities[1]->registerEventListener(this);
		}

		void CommunitySplitEvent::nextStep() {
			for (count com = 0; com < 2; ++com) {
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
						// FIXME: handle situation when not enough nodes could be sampled
						const std::vector<node> nodesToAdd = generator.communityNodeSampler.birthCommunityNodes(numNodesToAdd, communities[com]->getNodes());

						for (node u : nodesToAdd) {
							communities[com]->addNode(u);
						}
					}
				}
			}

			++currentStep;
			if (currentStep == numSteps) {
				active = false;
				communities[0]->setAvailable(true);
				communities[1]->setAvailable(true);
				communities[0]->unregisterEventListener(this);
				communities[1]->unregisterEventListener(this);
			}
		}

		void CommunitySplitEvent::notifyNodeRemovedFromCommunity(node u, CommunityPtr) {
			nodesToRemove[0].erase(u);
			nodesToRemove[1].erase(u);
		}

	}
}
