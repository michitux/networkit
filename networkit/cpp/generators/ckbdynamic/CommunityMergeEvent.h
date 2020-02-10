#ifndef COMMUNITY_MERGE_EVENT_H_
#define COMMUNITY_MERGE_EVENT_H_

#include "CommunityChangeEvent.h"

namespace NetworKit {
	namespace CKBDynamicImpl {
		class CKBDynamicImpl;

		/**
		 * Takes two communities.
		 * Add nodes from community A to B
		 * and vice-versa. However, while adding nodes,
		 * decrease the edge probability until it is so low
		 * that when combining both communities, we get the
		 * target edge probability (note: an existing overlap
		 * needs to be considered here). In the final step,
		 * just combine both edge sets and remove duplicate
		 * edges (needs special support in the community
		 * object).
		 */
		class CommunityMergeEvent : public CommunityChangeEvent {
		public:
			CommunityMergeEvent(CommunityPtr communityA, CommunityPtr communityB, count targetSize, count numSteps, CKBDynamicImpl& generator);
			virtual void nextStep() override;
			virtual void notifyNodeRemovedFromCommunity(node u, CommunityPtr com) override;
			virtual void notifyNodeAddedToCommunity(node u, CommunityPtr com) override;
			virtual bool canRemoveNode() const override;
		private:
			void mergeCommunities();
			std::array<Aux::SamplingSet<node>, 2> nodesToAddTo;
			std::array<CommunityPtr, 2> communities;
			count targetSize;
			// If the two communities were already merged into one.
			bool communitiesMerged;
		};
	}
}

#endif
