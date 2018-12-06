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

			double getProbability(count size) const;
		public:
			PowerlawCommunitySizeDistribution(count minSize, count maxSize, double gamma, double alpha, double densityGamma);

			virtual std::pair<count, double> drawCommunity() override;
		};
	}
}

#endif
