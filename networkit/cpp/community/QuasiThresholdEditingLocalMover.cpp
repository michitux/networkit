/*
 *
 */

#include <unordered_set>

#include <networkit/community/QuasiThresholdEditingLocalMover.hpp>
#include <networkit/generators/TreeReachabilityGraphGenerator.hpp>
#include <networkit/graph/DynamicForest.hpp>
#include <networkit/graph/GraphTools.hpp>
#include <networkit/auxiliary/SignalHandling.hpp>
#include <networkit/auxiliary/Log.hpp>
#include <networkit/auxiliary/Random.hpp>


NetworKit::QuasiThresholdEditingLocalMover::QuasiThresholdEditingLocalMover(
	const NetworKit::Graph &G, 
	const std::vector< NetworKit::node > &parent, 
	NetworKit::count maxIterations, 
	bool sortPaths,
	bool randomness,
	const std::vector<node> &order,
	count maxPlateauSize)
: G(G), 
	maxIterations(maxIterations), 
	sortPaths(sortPaths),
	randomness(randomness),
	order(order),
	maxPlateauSize(maxPlateauSize),
	numEdits(0),
	usedIterations(0),
	rootEqualBestParents(0) {
	forest = Graph(GraphTools::copyNodes(G), false, true);
	if(parent.size() == G.upperNodeIdBound()){
		//insert Run makes only sense if parents are trivial
		insertRun = 0;
		G.forNodes([&](node u) {
			if (parent[u] != none && parent[u] != u) {
				forest.addEdge(u, parent[u]);
			}
		});
	} else {
		insertRun = (order.size() == G.upperNodeIdBound());
	}
	runner = nullptr;
}

void NetworKit::QuasiThresholdEditingLocalMover::run() {
	/*
	Idea: for each node check if it can be inserted below any other node such that the number of edits is reduced.
	Additionally check if the node should become a root node and add any number of previous roots as children.
	For each child of the new parent check if it can be inserted below the moved node. Question: is this possible such that the total time for the move is linear?

	For a node that shall be moved, mark all neighbors in the original graph.
	For just checking if a node can be inserted as leaf under an anchor node or as parent of an anchor node only the neighbors of that anchor node need to be iterated
	and checked for the mark, then the number of necessary edits can be determined.

	Problem: In order to keep these checks easy we should remove the node from the graph. Therefore insertion in the previous position must remain possible. This means that
	we need to be able to add arbitrary children the node previously had. We find them when iterating over the neighbors of the parent.
	But can we determine which of them should be added?

	Another idea: Keep an efficient forest data structure (or simply rooted trees in a graph + a list of roots?) that allows walking up and down in the tree.
	Then do DFS in order to discover possible inserts/deletes.
	This should actually allow to determine which children to take. Furthermore then probably we don't need to store the inserted edges explicitly but we can determine
	them again whenever we delete a node as we iterate over its neighbors anyway and we can in the same turn also iterate over the tree and determine which neighbors were there
	and which were not by checking for marked nodes.
	*/
	
	runner = new EditingRunner(G, forest, maxIterations, sortPaths, randomness, order, maxPlateauSize, insertRun);
	runner->runLocalMover();
	forest = runner->getForest();
	usedIterations =  runner->getUsedIterations();
	numEdits = runner->getNumberOfEdits();
	plateauSize = runner->getPlateauSize();
	rootEqualBestParents = runner->getRootEqualBestParents();
	delete runner;
}



std::vector< NetworKit::node > NetworKit::QuasiThresholdEditingLocalMover::getParents() const {
	std::vector<node> parents(G.upperNodeIdBound());

	G.forNodes([&](node u) {
		node p = none;
		forest.forNeighborsOf(u, [&](node v) {
			p = v;
		});

		parents[u] = p;
	});

	return parents;
}

NetworKit::Graph NetworKit::QuasiThresholdEditingLocalMover::getQuasiThresholdGraph() const {
	TreeReachabilityGraphGenerator gen(forest);
	gen.run();
	return gen.getGraph();
}



NetworKit::count NetworKit::QuasiThresholdEditingLocalMover::getNumberOfEdits() const {
	return numEdits;
}




NetworKit::count NetworKit::QuasiThresholdEditingLocalMover::getUsedIterations() const {
	return usedIterations;
}

NetworKit::count NetworKit::QuasiThresholdEditingLocalMover::getPlateauSize() const {
	return plateauSize;
}

NetworKit::count NetworKit::QuasiThresholdEditingLocalMover::getRootEqualBestParents() const {
	return rootEqualBestParents;
}

