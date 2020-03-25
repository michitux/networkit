#ifndef QUASITHRESHOLDGREEDYBOUND_H
#define QUASITHRESHOLDGREEDYBOUND_H

#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {

class QuasiThresholdGreedyBound : public Algorithm {
public:
	QuasiThresholdGreedyBound(const Graph& G);
	void run() override;
	count getMinDistance() const { assureFinished(); return minDist; };
private:
	const Graph& G;
	count minDist;
};

} // namespace NetworKit

#endif // QUASITHRESHOLDGREEDYBOUND_H
