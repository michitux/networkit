#include "PowerlawCommunitySizeDistribution.h"
#include <tlx/unused.hpp>

namespace NetworKit {
	namespace CKBDynamicImpl {

		PowerlawCommunitySizeDistribution::PowerlawCommunitySizeDistribution(count minSize, count maxSize, double gamma, double alpha, double densityGamma, double minSplitRatio) : CommunitySizeDistribution(minSize, maxSize), sequence(minSize, maxSize, gamma), alpha(alpha), densityGamma(densityGamma), minSplitRatio(minSplitRatio) {
			sequence.run();
		}

		double PowerlawCommunitySizeDistribution::getProbability(count size) const {
			double prob = alpha / std::pow(size, densityGamma);
			if (prob > 1 || size <= 2) prob = 1;
			return prob;
		}

		std::pair<count, double> PowerlawCommunitySizeDistribution::drawCommunity() {
			count size = sequence.getDegree();
			double probability = getProbability(size);
			return std::make_pair(size, probability);
		}
	}

}
