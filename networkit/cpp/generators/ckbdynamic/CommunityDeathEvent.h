#ifndef COMMUNITY_DEATH_EVENT_H_
#define COMMUNITY_DEATH_EVENT_H_

#include "CommunityChangeEvent.h"

namespace NetworKit {
	namespace CKBDynamicImpl {
		class CKBDynamicImpl;

		/**
		 * Remove nodes over time.
		 * Lets community die at the end.
		 */
		class CommunityDeathEvent : public CommunityChangeEvent {
		public:
			CommunityDeathEvent(CommunityPtr community, count coreSize, count numSteps, CKBDynamicImpl& generator);
			virtual void nextStep() override;
		private:
			CommunityPtr community;
			count coreSize;
		};
	}
}

#endif
