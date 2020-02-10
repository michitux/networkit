#include "CommunitySplitEvent.h"

#include "CKBDynamicImpl.h"

namespace NetworKit {
	namespace CKBDynamicImpl {
		CommunitySplitEvent::CommunitySplitEvent(CommunityPtr community, count targetSizeA, count targetSizeB, count numSteps, CKBDynamicImpl& generator) : CommunityChangeEvent(generator, numSteps), targetSize({targetSizeA, targetSizeB}), communities({community, CommunityPtr(new Community(*community))}) {
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
			double fractionB = 1.0 - fractionA;
			std::array<count, 2> numNodesToRemove;
			numNodesToRemove[0] = fractionB * nodes.size();
			numNodesToRemove[1] = nodes.size() - numNodesToRemove[0];

			// We rounded a double - check that we did not round in the wrong direction.
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

			communities[0]->setCurrentEvent(this);
			communities[1]->setCurrentEvent(this);
		}

		void CommunitySplitEvent::nextStep() {
			assert(communities[0]->getNumberOfNodes() - nodesToRemove[0].size() <= targetSize[0]);
			assert(communities[1]->getNumberOfNodes() - nodesToRemove[1].size() <= targetSize[1]);

			for (count com = 0; com < 2; ++com) {
				assert(communities[com]->getNumberOfNodes() >= nodesToRemove[com].size());
				assert(communities[com]->getNumberOfNodes() - nodesToRemove[com].size() <= targetSize[com]);

				const count numNodesToRemove = std::ceil(static_cast<double>(nodesToRemove[com].size()) / (numSteps - currentStep));
				for (count i = 0; i < numNodesToRemove; ++i) {
					const node u = nodesToRemove[com].at(Aux::Random::index(nodesToRemove[com].size()));
					nodesToRemove[com].erase(u);

					assert(generator.hasNode(u));
					communities[com]->removeNode(u);
				}

				count numNodesToAdd = 0;

				// Add nodes to achieve target size, and at least the minimum size in the current step
				if (targetSize[com] > (communities[com]->getNumberOfNodes() - nodesToRemove[com].size())) {
					numNodesToAdd = std::ceil(static_cast<double>(targetSize[com] + nodesToRemove[com].size() - communities[com]->getNumberOfNodes()) / (numSteps - currentStep));
					if (communities[com]->getNumberOfNodes() + numNodesToAdd < generator.communitySizeSampler->getMinSize()) {
						numNodesToAdd = generator.communitySizeSampler->getMinSize() - communities[com]->getNumberOfNodes();
					}
				}

				communities[com]->setDesiredNumberOfNodes(communities[com]->getNumberOfNodes() + numNodesToAdd);
			}

			++currentStep;
			if (currentStep == numSteps) {
				active = false;
				for (count com = 0; com < 2; ++com) {
					assert(communities[com]->getDesiredNumberOfNodes() == targetSize[com]);
					communities[com]->setCurrentEvent(nullptr);
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

		bool CommunitySplitEvent::canRemoveNode() const {
			return false;
		}

	}
}
