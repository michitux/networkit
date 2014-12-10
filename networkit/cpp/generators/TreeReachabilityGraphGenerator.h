#ifndef TREEREACHABILITYGRAPHGENERATOR_H
#define TREEREACHABILITYGRAPHGENERATOR_H

#include "../graph/Graph.h"

namespace NetworKit {

class TreeReachabilityGraphGenerator {
public:
	TreeReachabilityGraphGenerator(const NetworKit::Graph &input);

	void run();

	Graph getGraph();
private:
	const Graph &input;
	Graph output;
};

} // namespace NetworKit

#endif // TREEREACHABILITYGRAPHGENERATOR_H
