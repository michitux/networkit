#ifndef COMMUNITY_BIRTH_EVENT_H_
#define COMMUNITY_BIRTH_EVENT_H_

#include "CommunityChangeEvent.h"

namespace NetworKit {
	namespace CKBDynamicImpl {
		class CKBDynamicImpl;

		/**
		 * Create community from a subset of given nodes.
		 * Adds further nodes over time.
		 */
		class CommunityBirthEvent : public CommunityChangeEvent {
		public:
			CommunityBirthEvent(count coreSize, count targetSize, count numSteps, CKBDynamicImpl& generator);
			CommunityBirthEvent(count numSteps, CKBDynamicImpl& generator);
			virtual void nextStep() override;
		private:
			count coreSize;
			count targetSize;
			CommunityPtr community;
		};
	}
}

#endif
