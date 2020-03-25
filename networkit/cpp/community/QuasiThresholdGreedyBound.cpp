/*
 *
 */

#include <networkit/community/QuasiThresholdGreedyBound.hpp>
#include <networkit/auxiliary/Log.hpp>
#include <robin_hood.h>

namespace NetworKit {

QuasiThresholdGreedyBound::QuasiThresholdGreedyBound(const Graph &G) : G(G), minDist(none) {

}

void QuasiThresholdGreedyBound::run() {
	struct NodePairHash {
		size_t operator()(const std::pair<node, node> &v) const {
			return robin_hood::hash_int((v.first<<32)+v.second);
		}
	};

	robin_hood::unordered_flat_set<std::pair<node, node>, NodePairHash> used_pairs;

	auto is_marked = [&used_pairs](node u, node v) {
		if (u > v) std::swap(u, v);
		return used_pairs.contains(std::make_pair(u, v));
	};

	auto mark = [&used_pairs,&is_marked](node u, node v) {
		assert(!is_marked(u, v));
		if (u > v) std::swap(u, v);
		used_pairs.emplace(u, v);
	};

	minDist = 0;
	std::vector<robin_hood::unordered_flat_set<node>> neighbors(G.upperNodeIdBound());

	G.forNodes([&](node u) {
		if (G.degree(u) < 2) return; // avoid degree-1-nodes as we skip them anyway below
		auto &uNeighbors = neighbors[u];
		uNeighbors.reserve(G.degree(u));
		G.forNeighborsOf(u, [&](node v) {
			uNeighbors.insert(v);
		});
	});

	// Find a neighbor of u that is not a neighbor of v and for which none of the node pairs is marked
	auto findUnmarkedNeighbor = [&](node u, node v) {
		for (node w : G.neighborRange(u)) {
			if (w != v && neighbors[v].count(w) == 0 && !is_marked(u, w) && !is_marked(v, w)) {
				return w;
			}
		}

		return none;
	};

	G.forEdges([&](node u, node v) {
		if (G.degree(u) < 2 || G.degree(v) < 2) return;
		if (is_marked(u, v)) return;

		// First check the lower-degree neighbor, then the higher-degree neighbor
		if (G.degree(u) > G.degree(v)) std::swap(u, v);

		node uNeighbor = findUnmarkedNeighbor(u, v);
		if (uNeighbor != none) {
			node vNeighbor = findUnmarkedNeighbor(v, u);

			if (vNeighbor != none) {
				mark(u, v);
				mark(u, uNeighbor);
				mark(v, uNeighbor);
				mark(u, vNeighbor);
				mark(v, vNeighbor);
				++minDist;
				if (minDist % 10000 == 0) {
					INFO("node ", u, " minDist increased to ", minDist);
				}
			}
		}
	});

	hasRun = true;
}

}
