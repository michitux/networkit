#include "CommunityDeathEvent.h"
#include "CKBDynamicImpl.h"

namespace NetworKit {
	namespace CKBDynamicImpl {

		CommunityDeathEvent::CommunityDeathEvent(CommunityPtr community, count coreSize, count numSteps, CKBDynamicImpl& generator) : CommunityChangeEvent(generator, numSteps), community(community), coreSize(coreSize) {
			community->setAvailable(false);
		}

		void CommunityDeathEvent::nextStep() {
			count numNodesToRemove;
			// Ensure that in the last step exactly coreSize nodes are removed.
			// If there are less than coreSize nodes remaining due to additional node
			// deletions let the community die earlier.
			if (currentStep < numSteps - 1 && community->getNumberOfNodes() >= coreSize) {
				numNodesToRemove = (community->getNumberOfNodes() - coreSize) / (numSteps - 1 - currentStep);
			} else {
				numNodesToRemove = community->getNumberOfNodes();
			}

			for (index i = 0; i < numNodesToRemove; ++i) {
				community->removeRandomNode();
			}

			community->setDesiredNumberOfNodes(community->getNumberOfNodes());

			++currentStep;

			if (currentStep == numSteps) {
				active = false;
				generator.removeCommunity(community);
			}
		}
	}
}
