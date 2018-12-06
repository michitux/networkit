#include "PowerlawCommunityMembershipDistribution.h"


namespace NetworKit {
	namespace CKBDynamicImpl {
		PowerlawCommunityMembershipDistribution::PowerlawCommunityMembershipDistribution(count minMemberships, count maxMemberships, double exponent) : sequence(minMemberships, maxMemberships, exponent) {
			sequence.run();
		}

		count PowerlawCommunityMembershipDistribution::drawMemberships() {
			return sequence.getDegree();
		}
	}
}
