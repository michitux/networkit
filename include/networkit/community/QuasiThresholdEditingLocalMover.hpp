#ifndef QUASITHRESHOLDEDITINGLOCALMOVER_H
#define QUASITHRESHOLDEDITINGLOCALMOVER_H

#include <map>
#include <networkit/community/QuasiThresholdMover/QuasiThresholdEditingLinear.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/graph/GraphTools.hpp>

namespace NetworKit {
namespace QuasiThresholdMoving {
class QuasiThresholdEditingLocalMover {
public:
    enum Initialization { TRIVIAL, EDITING, RANDOM_INSERT, ASC_DEGREE_INSERT, USER_DEFINED_INSERT };

    QuasiThresholdEditingLocalMover(const Graph &G, Initialization initializarion = TRIVIAL,
                                    count maxIterations = 5, bool sortPaths = true,
                                    bool randomness = false, count maxPlateauSize = 4,
                                    bool useBucketQueue = true);

    void run();

    Graph getQuasiThresholdGraph() const;
    count getNumberOfEdits() const;
    count getUsedIterations() const;
    count getPlateauSize() const;
    count getRootEqualBestParents() const;
    std::map<std::string, std::vector<count>> getRunningInfo() const;

    void setInsertionOrder(std::vector<node> order);

private:
    const Graph &G;
    Initialization initialization;
    count maxIterations;
    bool sortPaths;
    bool randomness;
    count maxPlateauSize;
    bool insertRun;
    bool useBucketQueue;

    std::vector<node> order;

    count usedIterations;
    count numEdits;
    count plateauSize;
    count rootEqualBestParents;

    Graph quasiThresholdGraph;

    std::map<std::string, std::vector<count>> runningInfo;
};
} // namespace QuasiThresholdMoving

} // namespace NetworKit

#endif // QUASITHRESHOLDEDITINGLOCALMOVER_H
