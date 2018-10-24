

#include "CommunityChangeEvent.h"

namespace NetworKit {
	namespace CKBDynamicImpl {
		CommunityChangeEvent::CommunityChangeEvent(CKBDynamicImpl& generator, count numSteps) : active(true), numSteps(numSteps), currentStep(0), generator(generator) {}

		void CommunityChangeEvent::adaptProbability(CommunityPtr com, double targetProb) {
			double prob = (com->getEdgeProbability() * (numSteps - currentStep - 1) + targetProb) / (numSteps - currentStep);
			com->changeEdgeProbability(prob);
		}

	}
}
