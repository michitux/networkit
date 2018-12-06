#ifndef CKBDYNAMIC_CUSTOM_COMMUNITY_MEMBERSHIP_DISTRIBUTION_H_
#define CKBDYNAMIC_CUSTOM_COMMUNITY_MEMBERSHIP_DISTRIBUTION_H_

#include "CommunityMembershipDistribution.h"

namespace NetworKit {
	class Graph;
	class Cover;
	namespace CKBDynamicImpl {
		class CustomCommunityMembershipDistribution : public CommunityMembershipDistribution {
		public:
			CustomCommunityMembershipDistribution(const Graph &G, const Cover &C);
			virtual count drawMemberships() override;
		private:
			std::vector<count> sequence;
		};
	}
}

#endif
