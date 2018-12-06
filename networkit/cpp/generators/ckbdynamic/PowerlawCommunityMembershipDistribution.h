#ifndef CKBDYNAMIC_POWERLAW_COMMUNITY_MEMBERSHIP_DISTRIBUTION_H_
#define CKBDYNAMIC_POWERLAW_COMMUNITY_MEMBERSHIP_DISTRIBUTION_H_

#include "CommunityMembershipDistribution.h"
#include "../PowerlawDegreeSequence.h"

namespace NetworKit {
	namespace CKBDynamicImpl {
		class PowerlawCommunityMembershipDistribution : public CommunityMembershipDistribution {
		public:
			PowerlawCommunityMembershipDistribution(count minMemberships, count maxMemberships, double exponent);
			virtual count drawMemberships() override;
		private:
			PowerlawDegreeSequence sequence;
		};
	}
}

#endif
