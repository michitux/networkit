#include "PowerlawCommunityMembershipDistribution.h"


namespace NetworKit {
	namespace CKBDynamicImpl {
		PowerlawCommunityMembershipDistribution::PowerlawCommunityMembershipDistribution(count minMemberships, count maxMemberships, double exponent) : CommunityMembershipDistribution(minMemberships, maxMemberships), sequence(minMemberships, maxMemberships, exponent) {
			sequence.run();
			avgMemberships = sequence.getExpectedAverageDegree();
		}

		count PowerlawCommunityMembershipDistribution::drawMemberships() {
			return sequence.getDegree();
		}
	}
}
