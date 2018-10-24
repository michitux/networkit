#ifndef COMMUNITY_CHANGE_EVENT_H_
#define COMMUNITY_CHANGE_EVENT_H_

#include "Community.h"

namespace NetworKit {
	namespace CKBDynamicImpl {
		class CKBDynamicImpl;

		class CommunityChangeEvent {
		public:
			CommunityChangeEvent(CKBDynamicImpl& generator, count numSteps);
			virtual void nextStep() = 0;
			bool isActive() { return active; };
		protected:
			void adaptProbability(CommunityPtr com, double targetProb);
			bool active;
			count numSteps;
			count currentStep;
			CKBDynamicImpl& generator;
		};
	}
}

#endif
