#include "CustomCommunitySizeDistribution.h"
#include "../../graph/Graph.h"
#include "../../structures/Cover.h"
#include <unordered_map>

namespace NetworKit {
	namespace CKBDynamicImpl {
		CustomCommunitySizeDistribution::CustomCommunitySizeDistribution(const Graph &G, const Cover &C) : CommunitySizeDistribution(std::numeric_limits<count>::max(), 0) {
			if (G.isWeighted()) throw std::runtime_error("Weighted graphs are not supported.");
			if (G.isDirected()) throw std::runtime_error("Directed graphs are not supported.");
			if (G.numberOfSelfLoops()) throw std::runtime_error("Graphs with self-loops are not supported.");

			std::unordered_map<index, std::pair<count, count>> communityNodesEdges;

			std::vector<index> shared;
			count globalEdges = 0;

			G.forNodes([&](node u) {
				std::set<index> comU(C.subsetsOf(u));

				for (index s : comU) {
					++communityNodesEdges[s].first;
				}

				G.forEdgesOf(u, [&](node v) {
					if (u > v) return;

					std::set<index> comV(C.subsetsOf(v));
					std::set_intersection(comU.begin(), comU.end(), comV.begin(), comV.end(), std::back_inserter(shared));

					if (shared.empty()) ++globalEdges;

					for (index s : shared) {
						++communityNodesEdges[s].second;
					}

					shared.clear();
				});
			});

			epsilon = globalEdges * 2.0 / (G.numberOfNodes() * (G.numberOfNodes() - 1));

			sequence.reserve(communityNodesEdges.size());
			for (auto it : communityNodesEdges) {
				const count size = it.second.first;
				const count edges = it.second.second;
				const double density = edges * 2.0 / (size * (size - 1));
				sequence.emplace_back(size, density);

				if (size < minSize) {
					minSize = size;
				}

				if (size > maxSize) {
					maxSize = size;
				}
			}

			if (sequence.empty()) throw std::runtime_error("Error, no communities found in reference.");
		}

		std::pair<count, double> CustomCommunitySizeDistribution::drawCommunity() {
			return Aux::Random::choice(sequence);
		}

		double CustomCommunitySizeDistribution::getEpsilon() {
			return epsilon;
		}

	}
}
