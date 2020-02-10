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

		public:
			PowerlawCommunitySizeDistribution(count minSize, count maxSize, double gamma, double alpha, double densityGamma);

			virtual count drawCommunitySize() override;
			virtual double getCommunityDensity(count size) override;
		};
	}
}

#endif
