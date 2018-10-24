#include "CommunityBirthEvent.h"
#include "CKBDynamicImpl.h"

namespace NetworKit {
	namespace CKBDynamicImpl {

	CommunityBirthEvent::CommunityBirthEvent(count coreSize, count targetSize, double edgeProbability, count numSteps, CKBDynamicImpl& generator) : CommunityChangeEvent(generator, numSteps), coreSize(coreSize), targetSize(targetSize), community(new Community(edgeProbability, generator)) {}

	CommunityBirthEvent::CommunityBirthEvent(count numSteps, CKBDynamicImpl& generator) : CommunityChangeEvent(generator, numSteps)	{
		double edgeProbability;
		std::tie(targetSize, edgeProbability) = generator.communitySizeSampler->drawCommunity();
		coreSize = std::max<count>(0.1 * targetSize, generator.communitySizeSampler->getMinSize());
		community = CommunityPtr(new Community(edgeProbability, generator));
	}

	void CommunityBirthEvent::nextStep() {
		count numNodesToAdd = (targetSize - community->getNumberOfNodes()) / (numSteps - currentStep);

		// Ensure that after every step the community has at least coreSize nodes,
		// even if nodes have been deleted in the meantime.
		if (community->getNumberOfNodes() < coreSize) {
			numNodesToAdd = std::max(coreSize - community->getNumberOfNodes(), numNodesToAdd);
		}

		std::vector<node> nodesToAdd = generator.communityNodeSampler.birthCommunityNodes(numNodesToAdd, community->getNodes());

		// If we did not get enough nodes to make the community larger than the minimum size
		// let the community die.
		if (nodesToAdd.size() + community->getNumberOfNodes() < coreSize) {
			INFO("Sampler returned too few nodes, letting the community die.");
			active = false;
			while (community->getNumberOfNodes() > 0) {
				community->removeRandomNode();
			}
			generator.removeCommunity(community);
		} else {
			for (node u : nodesToAdd) {
				assert(generator.hasNode(u));
				community->addNode(u);
			}

			++currentStep;

			if (currentStep == numSteps) {
				if (community->getNumberOfNodes() == targetSize) {
					active = false;
					community->setAvailable(true);
				} else {
					// Delay event completion because not enough nodes could be requested.
					INFO("Delaying community birth because not enough nodes could be sampled.");
					--currentStep;
				}
			}
		}
	}
	}
}
