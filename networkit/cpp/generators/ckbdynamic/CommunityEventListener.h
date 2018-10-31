#ifndef CKBDYNAMIC_COMMUNITY_EVENT_LISTENER_H_
#define CKBDYNAMIC_COMMUNITY_EVENT_LISTENER_H_

#include "Community.h"

namespace NetworKit {
	namespace CKBDynamicImpl {
		class CommunityEventListener {
		public:
			virtual ~CommunityEventListener() = default;
			virtual void notifyNodeAddedToCommunity(node, CommunityPtr);
			virtual void notifyNodeRemovedFromCommunity(node, CommunityPtr);
			virtual void notifyEdgeAddedToCommunity(node, node, CommunityPtr);
			virtual void notifyEdgeRemovedFromCommunity(node, node, CommunityPtr);
		};
	};
}
#endif
