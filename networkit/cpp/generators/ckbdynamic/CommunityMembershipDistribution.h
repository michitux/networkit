#ifndef CKBDYNAMIC_COMMUNITY_MEMBERSHIP_DISTRIBUTION_H_
#define CKBDYNAMIC_COMMUNITY_MEMBERSHIP_DISTRIBUTION_H_

#include "../../Globals.h"

namespace NetworKit {
	namespace CKBDynamicImpl {
		class CommunityMembershipDistribution {
		public:
			virtual count drawMemberships() = 0;
		};
	}
}

#endif
