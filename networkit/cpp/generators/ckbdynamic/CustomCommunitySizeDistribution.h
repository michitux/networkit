#ifndef CKBDYNAMIC_CUSTOM_COMMUNITY_SIZE_DISTRIBUTION_H_
#define CKBDYNAMIC_CUSTOM_COMMUNITY_SIZE_DISTRIBUTION_H_

#include "CommunitySizeDistribution.h"

namespace NetworKit {
	class Graph;
	class Cover;
	namespace CKBDynamicImpl {
		class CustomCommunitySizeDistribution : public CommunitySizeDistribution {
		public:
			CustomCommunitySizeDistribution(const Graph &G, const Cover &C);

			virtual count drawCommunitySize() override;
			virtual double getCommunityDensity(count size) override;
			double getEpsilon();
		private:
			std::vector<std::pair<count, double>> sequence;
			double epsilon;
		};
	}
}

#endif
