
#include <networkit/community/QuasiThresholdEditingLocalMover.hpp>
#include <networkit/community/QuasiThresholdMover/EditingRunner.hpp>

namespace NetworKit {
namespace QuasiThresholdMoving {
QuasiThresholdEditingLocalMover::QuasiThresholdEditingLocalMover(
    const Graph &G, Initialization initialization, count maxIterations, bool sortPaths,
    bool randomness, count maxPlateauSize, bool useBucketQueue)
    : G(G), initialization(initialization), maxIterations(maxIterations), sortPaths(sortPaths),
      randomness(randomness), maxPlateauSize(maxPlateauSize), useBucketQueue(useBucketQueue),
      numEdits(0), usedIterations(0), rootEqualBestParents(0) {}

void QuasiThresholdEditingLocalMover::run() {
    EditingRunner runner(G, initialization, maxIterations, sortPaths, randomness, maxPlateauSize,
                         useBucketQueue, order);
    runner.runLocalMover();
    usedIterations = runner.getUsedIterations();
    numEdits = runner.getNumberOfEdits();
    plateauSize = runner.getPlateauSize();
    rootEqualBestParents = runner.getRootEqualBestParents();
    quasiThresholdGraph = runner.getQuasiThresholdGraph();
    runningInfo = runner.getRunningInfo();
}

Graph QuasiThresholdEditingLocalMover::getQuasiThresholdGraph() const {
    return quasiThresholdGraph;
}

count QuasiThresholdEditingLocalMover::getNumberOfEdits() const {
    return numEdits;
}

count QuasiThresholdEditingLocalMover::getUsedIterations() const {
    return usedIterations;
}

count QuasiThresholdEditingLocalMover::getPlateauSize() const {
    return plateauSize;
}

count QuasiThresholdEditingLocalMover::getRootEqualBestParents() const {
    return rootEqualBestParents;
}

std::map<std::string, std::vector<count>> QuasiThresholdEditingLocalMover::getRunningInfo() const {
    return runningInfo;
}

void QuasiThresholdEditingLocalMover::setInsertionOrder(std::vector<node> order) {
    this->order = order;
}

} // namespace QuasiThresholdMoving
} // namespace NetworKit
