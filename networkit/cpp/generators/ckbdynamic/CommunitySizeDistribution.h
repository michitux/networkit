#ifndef CKBDYNAMIC_COMMUNITY_SIZE_DISTRIBUTION_H_
#define CKBDYNAMIC_COMMUNITY_SIZE_DISTRIBUTION_H_

#include "../../Globals.h"

namespace NetworKit {
	namespace CKBDynamicImpl {
		class CommunitySizeDistribution {
		public:
			CommunitySizeDistribution(count minSize, count maxSize) : minSize(minSize), maxSize(maxSize) {};
			virtual std::pair<count, double> drawCommunity() = 0;
			virtual std::pair<count, double> mergeCommunities(count sizeA, double probabilityA, count sizeB, double probabilityB, count combinedNodes) = 0;
			virtual std::pair<std::pair<count, double>, std::pair<count, double>> splitCommunity(count size, double probability) = 0;
			count getMinSize() const { return minSize; };
			count getMaxSize() const { return maxSize; };
		protected:
			count minSize;
			count maxSize;
		};
	}
}

#endif
