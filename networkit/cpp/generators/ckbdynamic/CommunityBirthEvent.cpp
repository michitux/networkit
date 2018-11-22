#include "CommunityBirthEvent.h"
#include "CKBDynamicImpl.h"

namespace NetworKit {
	namespace CKBDynamicImpl {

	CommunityBirthEvent::CommunityBirthEvent(count coreSize, count targetSize, double edgeProbability, count numSteps, CKBDynamicImpl& generator) : CommunityChangeEvent(generator, numSteps), coreSize(coreSize), targetSize(targetSize), community(new Community(edgeProbability, generator)) {
		community->setCurrentEvent(this);
	}

	CommunityBirthEvent::CommunityBirthEvent(count numSteps, CKBDynamicImpl& generator) : CommunityChangeEvent(generator, numSteps)	{
		double edgeProbability;
		std::tie(targetSize, edgeProbability) = generator.communitySizeSampler->drawCommunity();
		coreSize = std::max<count>(0.1 * targetSize, generator.communitySizeSampler->getMinSize());
		community = CommunityPtr(new Community(edgeProbability, generator));
		community->setCurrentEvent(this);
	}

	void CommunityBirthEvent::nextStep() {
		count nextSize = targetSize / (numSteps - currentStep);

		// Ensure that after every step the community has at least coreSize nodes,
		// even if nodes have been deleted in the meantime.
		if (nextSize < coreSize) {
			nextSize = coreSize;
		}

		community->setDesiredNumberOfNodes(nextSize);
		// If we did not get enough nodes to make the community larger than the minimum size
		// let the community die.
		++currentStep;

		if (currentStep == numSteps) {
			active = false;
			community->setCurrentEvent(nullptr);
		}
	}
	}
}
