#ifndef COMMUNITY_CHANGE_EVENT_H_
#define COMMUNITY_CHANGE_EVENT_H_

#include "Community.h"

namespace NetworKit {
	namespace CKBDynamicImpl {
		class CKBDynamicImpl;

		class CommunityChangeEvent {
		public:
			CommunityChangeEvent(CKBDynamicImpl& generator, count numSteps);
			virtual ~CommunityChangeEvent() = default;
			virtual void nextStep() = 0;
			virtual void notifyNodeAddedToCommunity(node, CommunityPtr);
			virtual void notifyNodeRemovedFromCommunity(node, CommunityPtr);
			bool isActive() { return active; };
			virtual bool canRemoveNode() const;
		protected:
			bool active;
			count numSteps;
			count currentStep;
			CKBDynamicImpl& generator;
		};
	}
}

#endif
