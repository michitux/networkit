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
			CommunitySplitEvent(CommunityPtr community, count targetSizeA, double targetEdgeProbabilityA, count targetSizeB, double targetEdgeProbabilityB, count numSteps, CKBDynamicImpl& generator);
			virtual void nextStep() override;
		private:
			std::array<std::vector<node>, 2> nodesToRemove;
			std::array<double, 2> targetEdgeProbability;
			std::array<CommunityPtr, 2> communities;
		};
	}
}

#endif
