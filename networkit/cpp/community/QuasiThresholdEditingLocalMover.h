#ifndef QUASITHRESHOLDEDITINGLOCALMOVER_H
#define QUASITHRESHOLDEDITINGLOCALMOVER_H

#include "../graph/Graph.h"

namespace NetworKit {

class QuasiThresholdEditingLocalMover {
public:
	QuasiThresholdEditingLocalMover(const NetworKit::Graph &G, const std::vector< NetworKit::node > &parent, NetworKit::count maxIterations);

	void run();

	Graph getQuasiThresholdGraph();
	count getNumberOfEdits();
private:
	const Graph& G;
	Graph forest;
	count maxIterations;
	count numEdits;
};

} // namespace NetworKit

#endif // QUASITHRESHOLDEDITINGLOCALMOVER_H
