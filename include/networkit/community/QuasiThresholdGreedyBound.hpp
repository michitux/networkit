#ifndef QUASITHRESHOLDGREEDYBOUND_H
#define QUASITHRESHOLDGREEDYBOUND_H

#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {

/*
 * Simple greedy bound for the quasi-threshold editing distance.
 *
 * Computes a bound by greedily packing forbidden subgraphs as described in
 *
 * Gottesbüren, L., Hamann, M., Schoch, P., Strasser, B., Wagner, D., &
 * Zühlsdorf, S. (2020). Engineering Exact Quasi-Threshold Editing. In S. Faro
 * & D. Cantone (Eds.), 18th International Symposium on Experimental Algorithms
 * (SEA 2020) (Vol. 160, p. 10:1–10:14). Schloss Dagstuhl–Leibniz-Zentrum für
 * Informatik. https://doi.org/10.4230/LIPIcs.SEA.2020.10
 */
class QuasiThresholdGreedyBound : public Algorithm {
public:
	/**
	 * Initialize the greedy quasi-threshold editing bound
	 *
	 * @param G The graph to compute the bound for
	 */
	QuasiThresholdGreedyBound(const Graph& G);
	void run() override;
	/**
	 * Get the minimum quasi-threshold editing distance
	 *
	 * @return The minimum editing distance
	 */
	count getMinDistance() const { assureFinished(); return minDist; };
private:
	const Graph& G;
	count minDist;
};

} // namespace NetworKit

#endif // QUASITHRESHOLDGREEDYBOUND_H
