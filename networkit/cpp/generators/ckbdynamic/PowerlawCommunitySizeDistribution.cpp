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

		std::pair<count, double> PowerlawCommunitySizeDistribution::mergeCommunities(count sizeA, double probabilityA, count sizeB, double probabilityB, count combinedNodes) {
			tlx::unused(sizeA);
			tlx::unused(probabilityA);
			tlx::unused(sizeB);
			tlx::unused(probabilityB);
			//assert(getProbability(sizeA) == probabilityA);
			//assert(getProbability(sizeB) == probabilityB);

			return std::make_pair(combinedNodes, getProbability(combinedNodes));
		}

		std::pair<std::pair<count, double>, std::pair<count, double>> PowerlawCommunitySizeDistribution::splitCommunity(count size, double) {
			count minSize = std::max<count>(minSplitRatio * size, sequence.getMinimumDegree());
			if (minSize * 2 > size) {
				throw std::runtime_error("Given community too small to split");
			}

			std::mt19937 gen(rand());
			std::uniform_int_distribution<> dis(minSize, size - minSize);
			int sizeA = dis(gen);
			int sizeB = size-sizeA;

			return std::make_pair(std::make_pair(sizeA, getProbability(sizeA)),std::make_pair(sizeB, getProbability(sizeB)));
		}
	}

}
