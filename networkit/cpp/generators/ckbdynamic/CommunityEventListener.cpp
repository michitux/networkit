#include "CommunityEventListener.h"


namespace NetworKit {
	namespace CKBDynamicImpl {
		void CommunityEventListener::notifyNodeAddedToCommunity(node, CommunityPtr) {}
		void CommunityEventListener::notifyNodeRemovedFromCommunity(node, CommunityPtr) {}
		void CommunityEventListener::notifyEdgeAddedToCommunity(node, node, CommunityPtr) {}
		void CommunityEventListener::notifyEdgeRemovedFromCommunity(node, node, CommunityPtr) {}
	}
}
