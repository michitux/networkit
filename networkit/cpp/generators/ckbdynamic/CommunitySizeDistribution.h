#ifndef CKBDYNAMIC_COMMUNITY_SIZE_DISTRIBUTION_H_
#define CKBDYNAMIC_COMMUNITY_SIZE_DISTRIBUTION_H_

#include "../../Globals.h"

namespace NetworKit {
	namespace CKBDynamicImpl {
		class CommunitySizeDistribution {
		public:
			CommunitySizeDistribution(count minSize, count maxSize) : minSize(minSize), maxSize(maxSize) {};
			virtual ~CommunitySizeDistribution() = default;
			virtual count drawCommunitySize() = 0;
			virtual double getCommunityDensity(count size) = 0;
			count getMinSize() const { return minSize; };
			double getAverageSize() const { return avgSize; };
			count getMaxSize() const { return maxSize; };
		protected:
			count minSize;
			double avgSize;
			count maxSize;
		};
	}
}

#endif
