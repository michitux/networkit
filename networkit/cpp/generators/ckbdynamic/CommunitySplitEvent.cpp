#include "CommunitySplitEvent.h"

#include "CKBDynamicImpl.h"

namespace NetworKit {
	namespace CKBDynamicImpl {
		CommunitySplitEvent::CommunitySplitEvent(CommunityPtr community, count targetSizeA, double targetEdgeProbabilityA, count targetSizeB, double targetEdgeProbabilityB, count numSteps, CKBDynamicImpl& generator) : CommunityChangeEvent(generator, numSteps), targetEdgeProbability({targetEdgeProbabilityA, targetEdgeProbabilityB}), communities({community, CommunityPtr(new Community(*community))}) {
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
				nodesToRemove[0].push_back(u);
				nodesToRemove[1].push_back(u);
			}

			// Calculate how many nodes need to be removed from the two communities.
			// Due to the previous loop, both numbers can sum up to at most nodes.size().
			// If they sum up to less, we just keep some nodes in the overlap.
			const std::array<count, 2> numNodesToRemove {nodes.size() - targetSizeA, nodes.size() - targetSizeB};

			for (count com = 0; com < 2; ++com) {
				for (count i = 0; i < numNodesToRemove[com]; ++i) {
					node u = nodes.back();
					nodes.pop_back();
					nodesToRemove[com].push_back(u);
				}
			}

			communities[0]->setAvailable(false);
			communities[1]->setAvailable(false);
		}

		void CommunitySplitEvent::nextStep() {
			for (count com = 0; com < 2; ++com) {
				const count numNodesToRemove = nodesToRemove[com].size() / (numSteps - currentStep);
				for (count i = 0; i < numNodesToRemove; ++i) {
					const node u = nodesToRemove[com].back();
					nodesToRemove[com].pop_back();
					if (generator.hasNode(u)) {
						communities[com]->removeNode(u);
					}
				}

				adaptProbability(communities[com], targetEdgeProbability[com]);
			}

			++currentStep;
			if (currentStep == numSteps) {
				active = false;
				communities[0]->setAvailable(true);
				communities[1]->setAvailable(true);
			}
		}

	}
}
