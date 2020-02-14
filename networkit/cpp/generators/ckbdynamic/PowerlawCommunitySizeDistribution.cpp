#include "PowerlawCommunitySizeDistribution.h"
#include <tlx/unused.hpp>

namespace NetworKit {
	namespace CKBDynamicImpl {

		PowerlawCommunitySizeDistribution::PowerlawCommunitySizeDistribution(count minSize, count maxSize, double gamma, double alpha, double densityGamma) : CommunitySizeDistribution(minSize, maxSize), sequence(minSize, maxSize, gamma), alpha(alpha), densityGamma(densityGamma) {
			sequence.run();
			avgSize = sequence.getExpectedAverageDegree();
		}

		double PowerlawCommunitySizeDistribution::getCommunityDensity(count size) {
			double prob = alpha / std::pow(size, densityGamma);
			if (prob > 1 || size <= 2) prob = 1;
			return prob;
		}

		count PowerlawCommunitySizeDistribution::drawCommunitySize() {
			return sequence.getDegree();
		}
	}

}
