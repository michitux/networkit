/*
 *
 */

#include <unordered_set>

#include <networkit/community/QuasiThresholdEditingLocalMover.hpp>
#include <networkit/community/DynamicForest.hpp>
#include <networkit/graph/GraphTools.hpp>
#include <networkit/auxiliary/SignalHandling.hpp>
#include <networkit/auxiliary/Log.hpp>
#include <networkit/auxiliary/Random.hpp>


NetworKit::QuasiThresholdEditingLocalMover::QuasiThresholdEditingLocalMover(
	const NetworKit::Graph &G, 
	Initialization initialization, 
	NetworKit::count maxIterations, 
	bool sortPaths,
	bool randomness,
	count maxPlateauSize,
	bool useBucketQueue)
: G(G), 
	initialization(initialization),
	maxIterations(maxIterations), 
	sortPaths(sortPaths),
	randomness(randomness),
	maxPlateauSize(maxPlateauSize),
	useBucketQueue(useBucketQueue),
	numEdits(0),
	usedIterations(0),
	rootEqualBestParents(0) {
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
	runner = new EditingRunner(G, initialization, maxIterations, sortPaths, randomness, maxPlateauSize, useBucketQueue, order);
	runner->runLocalMover();
	usedIterations =  runner->getUsedIterations();
	numEdits = runner->getNumberOfEdits();
	plateauSize = runner->getPlateauSize();
	rootEqualBestParents = runner->getRootEqualBestParents();
	quasiThresholdGraph = runner->getQuasiThresholdGraph();
	delete runner;
}


NetworKit::Graph NetworKit::QuasiThresholdEditingLocalMover::getQuasiThresholdGraph() const {
	return quasiThresholdGraph;
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

void NetworKit::QuasiThresholdEditingLocalMover::setInsertionOrder(std::vector<node> order) {
	this->order = order;
}


