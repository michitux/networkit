#include "CustomCommunityMembershipDistribution.h"
#include "../../graph/Graph.h"
#include "../../structures/Cover.h"

namespace NetworKit {
	namespace CKBDynamicImpl {
		CustomCommunityMembershipDistribution::CustomCommunityMembershipDistribution(const Graph &G, const Cover &C) {
			sequence.reserve(G.numberOfNodes());
			G.forNodes([&](node u) {
					   sequence.push_back(C.subsetsOf(u).size());
				   });
		}

		count CustomCommunityMembershipDistribution::drawMemberships() {
			return Aux::Random::choice(sequence);
		}
	}
}
