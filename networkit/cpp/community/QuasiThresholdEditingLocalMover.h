#ifndef QUASITHRESHOLDEDITINGLOCALMOVER_H
#define QUASITHRESHOLDEDITINGLOCALMOVER_H

#include "../graph/Graph.h"

namespace NetworKit {

class QuasiThresholdEditingLocalMover {
public:
	QuasiThresholdEditingLocalMover(const NetworKit::Graph &G, const std::vector< NetworKit::node > &parent, NetworKit::count maxIterations, bool moveSubtrees = false);

	void run();

	Graph getQuasiThresholdGraph();
	count getNumberOfEdits();
	count getUsedIterations() const;
	std::vector<node> getParents() const;
private:
	count countNumberOfEdits() const;
	const Graph& G;
	Graph forest;
	count maxIterations;
	count usedIterations;
	bool moveSubtrees;
	count numEdits;
};

} // namespace NetworKit

#endif // QUASITHRESHOLDEDITINGLOCALMOVER_H
