#ifndef QUASITHRESHOLDEDITINGLOCALMOVER_H
#define QUASITHRESHOLDEDITINGLOCALMOVER_H

#include <networkit/graph/Graph.hpp>
#include <networkit/structures/Cover.hpp>

namespace NetworKit {

class QuasiThresholdEditingLocalMover {
public:
	QuasiThresholdEditingLocalMover(const NetworKit::Graph &G, const std::vector< NetworKit::node > &parent, NetworKit::count maxIterations, bool moveSubtrees = false);

	void run();

	Graph getQuasiThresholdGraph() const;
	count getNumberOfEdits() const;
	count getUsedIterations() const;
	std::vector<node> getParents() const;
	Cover getCover(NetworKit::index mergeDepth) const;
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
