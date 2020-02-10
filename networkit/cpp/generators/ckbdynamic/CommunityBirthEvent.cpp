#include "CommunityBirthEvent.h"
#include "CKBDynamicImpl.h"

namespace NetworKit {
	namespace CKBDynamicImpl {

	CommunityBirthEvent::CommunityBirthEvent(count coreSize, count targetSize, count numSteps, CKBDynamicImpl& generator) : CommunityChangeEvent(generator, numSteps), coreSize(coreSize), targetSize(targetSize), community(new Community(generator)) {
		community->setCurrentEvent(this);
	}

	CommunityBirthEvent::CommunityBirthEvent(count numSteps, CKBDynamicImpl& generator) : CommunityChangeEvent(generator, numSteps)	{
		targetSize = generator.communitySizeSampler->drawCommunitySize();
		coreSize = generator.communitySizeSampler->getMinSize();
		community = CommunityPtr(new Community(generator));
		community->setCurrentEvent(this);
	}

	void CommunityBirthEvent::nextStep() {
		count nextSize;
		if (currentStep + 1 == numSteps) {
			nextSize = targetSize;
		} else if (currentStep == 0) {
			nextSize = coreSize;
		} else {
			nextSize = community->getDesiredNumberOfNodes() +
				std::ceil(static_cast<double>(targetSize - community->getDesiredNumberOfNodes()) / (numSteps - currentStep));
			assert(nextSize <= targetSize);
                }

		community->setDesiredNumberOfNodes(nextSize);

		++currentStep;

		if (currentStep == numSteps || nextSize == targetSize) {
			active = false;
			community->setCurrentEvent(nullptr);
		}
	}
	}
}
