#ifndef QUASITHRESHOLDEDITINGLINEAR_H
#define QUASITHRESHOLDEDITINGLINEAR_H

#include "../graph/Graph.h"

namespace NetworKit {

class QuasiThresholdEditingLinear {
public:
	QuasiThresholdEditingLinear(const Graph& G);

	void run();

	std::vector<node> getParents() const;
	Graph getDefiningForest() const;
	Graph getQuasiThresholdGraph() const;

private:
	const Graph& G;
	std::vector<node> parent;
	bool hasRun;
};

} // namespace NetworKit

#endif // QUASITHRESHOLDEDITINGLINEAR_H
