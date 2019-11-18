/*
 *
 */

#include <stack>

#include <networkit/generators/TreeReachabilityGraphGenerator.hpp>

NetworKit::TreeReachabilityGraphGenerator::TreeReachabilityGraphGenerator(const NetworKit::Graph &input) : input(input) {

}


void NetworKit::TreeReachabilityGraphGenerator::run() {
	output = Graph(input.upperNodeIdBound());

	// if necessary, delete deleted nodes from the input.
	if (input.upperNodeIdBound() > input.numberOfNodes()) {
		for (node u = 0; u < input.upperNodeIdBound(); ++u) {
			if (!input.hasNode(u)) {
				output.removeNode(u);
			}
		}
	}

	input.forNodes([&](node u) {
		node curNode = u;

		while (curNode != none) {
			node nextNode = none;
			input.forNeighborsOf(curNode, [&](node v) {
				output.addEdge(u, v);
				nextNode = v;
			});
			curNode = nextNode;
		}
	});
}

NetworKit::Graph NetworKit::TreeReachabilityGraphGenerator::getGraph() {
	return output;
}
