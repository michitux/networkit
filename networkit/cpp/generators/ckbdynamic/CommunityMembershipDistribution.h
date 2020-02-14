#ifndef CKBDYNAMIC_COMMUNITY_MEMBERSHIP_DISTRIBUTION_H_
#define CKBDYNAMIC_COMMUNITY_MEMBERSHIP_DISTRIBUTION_H_

#include "../../Globals.h"

namespace NetworKit {
	namespace CKBDynamicImpl {
		class CommunityMembershipDistribution {
		public:
			CommunityMembershipDistribution(count minMemberships, count maxMemberships) : minMemberships(minMemberships), maxMemberships(maxMemberships) {};
			virtual ~CommunityMembershipDistribution() = default;
			virtual count drawMemberships() = 0;
			count getMinimumMemberships() const { return minMemberships; }
			double getAverageMemberships() const { return avgMemberships; };
			count getMaximumMemberships() const { return maxMemberships; }
		protected:
			count minMemberships;
			double avgMemberships;
			count maxMemberships;
		};
	}
}

#endif
