#ifndef CKBDYNAMIC_POWERLAW_COMMUNITY_SIZE_DISTRIBUTION_H_
#define CKBDYNAMIC_POWERLAW_COMMUNITY_SIZE_DISTRIBUTION_H_

#include "../../Globals.h"
#include "CommunitySizeDistribution.h"
#include "../PowerlawDegreeSequence.h"

namespace NetworKit {
	namespace CKBDynamicImpl {
		class PowerlawCommunitySizeDistribution : public CommunitySizeDistribution {
		private:
			PowerlawDegreeSequence sequence;
			double alpha;
			double densityGamma;
			double minSplitRatio;

			double getProbability(count size) const;
		public:
			PowerlawCommunitySizeDistribution(count minSize, count maxSize, double gamma, double alpha, double densityGamma, double minSplitRatio);

			virtual std::pair<count, double> drawCommunity() override;

			virtual std::pair<count, double> mergeCommunities(count sizeA, double probabilityA, count sizeB, double probabilityB, count combinedNodes) override;

			virtual std::pair<std::pair<count, double>, std::pair<count, double>> splitCommunity(count size, double probability) override;
		};
	}
}

#endif
