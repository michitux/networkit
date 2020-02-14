#include "CustomCommunityMembershipDistribution.h"
#include "../../graph/Graph.h"
#include "../../structures/Cover.h"

namespace NetworKit {
	namespace CKBDynamicImpl {
		CustomCommunityMembershipDistribution::CustomCommunityMembershipDistribution(const Graph &G, const Cover &C) : CommunityMembershipDistribution(std::numeric_limits<count>::max(), 0) {
			sequence.reserve(G.numberOfNodes());
			count sum = 0;
			G.forNodes([&](node u) {
				count numSubsets = C.subsetsOf(u).size();
				if (numSubsets < minMemberships) {
					minMemberships = numSubsets;
				}

				if (numSubsets > maxMemberships) {
					maxMemberships = numSubsets;
				}

				sum += numSubsets;
				sequence.push_back(numSubsets);
			});

			avgMemberships = sum * 1.0 / G.numberOfNodes();
		}

		count CustomCommunityMembershipDistribution::drawMemberships() {
			return Aux::Random::choice(sequence);
		}
	}
}
