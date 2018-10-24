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
			CommunityMergeEvent(CommunityPtr communityA, CommunityPtr communityB, count numSteps, CKBDynamicImpl& generator);
			virtual void nextStep() override;
		private:
			std::vector<node> nodesToAddToA;
			std::vector<node> nodesToAddToB;
			double targetEdgeProbability;
			double targetEdgeProbabilityPerCommunity;
			CommunityPtr communityA;
			CommunityPtr communityB;
		};
	}
}

#endif
