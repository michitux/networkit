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
	std::vector<bool> uMarker(G.upperNodeIdBound(), false);

	G.forEdges([&](node u, node v) {
		if (G.degree(u) < 2 || G.degree(v) < 2) return;
		if (is_marked(u, v)) return;

		if (G.degree(u) > G.degree(v)) std::swap(u, v);

		G.forNeighborsOf(u, [&](node w) {
			uMarker[w] = (w != v);
		});

		node vNeighbor = none;

		G.forNeighborsOf(v, [&](node w) {
			if (uMarker[w]) {
				uMarker[w] = false;
			} else if (w != u && vNeighbor == none && !(is_marked(u, w) || is_marked(v, w))) {
				vNeighbor = w;
			}
		});

		if (vNeighbor != none) {
			bool found = false;
			G.forNeighborsOf(u, [&](node w) {
				if (found) return;

				if (uMarker[w] && (!(is_marked(u, w) || is_marked(v, w)))) {
					mark(u, v);
					mark(u, w);
					mark(v, w);
					mark(u, vNeighbor);
					mark(v, vNeighbor);
					found = true;
					++minDist;
					if (minDist % 1024 == 0) {
						INFO("node ", u, " minDist increased to ", minDist);
					}
				}
			});
		}

		G.forNeighborsOf(u, [&](node w) {
			uMarker[w] = false;
		});
	});

	hasRun = true;
}

}
