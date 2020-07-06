
#include <networkit/community/QuasiThresholdEditingLocalMover.hpp>


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

