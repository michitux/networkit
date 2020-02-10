#ifndef COMMUNITY_SPLIT_EVENT_H_
#define COMMUNITY_SPLIT_EVENT_H_

#include "CommunityChangeEvent.h"

namespace NetworKit {
	namespace CKBDynamicImpl {
		class CKBDynamicImpl;

		/**
		 * Takes one community.
		 * Creates a second community as exact copy with the same edge set.
		 * Removes parts of both communities over time.
		 * Increases the density of both communities over time.
		 */
		class CommunitySplitEvent : public CommunityChangeEvent {
		public:
			CommunitySplitEvent(CommunityPtr community, count targetSizeA, count targetSizeB, count numSteps, CKBDynamicImpl& generator);
			virtual void nextStep() override;
			virtual void notifyNodeRemovedFromCommunity(node u, CommunityPtr com) override;
			virtual bool canRemoveNode() const override;
		private:
			std::array<Aux::SamplingSet<node>, 2> nodesToRemove;
                        const std::array<count, 2> targetSize;
                        const std::array<CommunityPtr, 2> communities;
		};
	}
}

#endif
