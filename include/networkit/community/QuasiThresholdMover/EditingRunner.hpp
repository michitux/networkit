#ifndef EDITINGRUNNER_H
#define EDITINGRUNNER_H

#include <limits>
#include <map>
#include <networkit/auxiliary/PerfEventCountHardware.hpp>
#include <networkit/auxiliary/SignalHandling.hpp>
#include <networkit/auxiliary/Timer.hpp>
#include <networkit/community/QuasiThresholdEditingLocalMover.hpp>
#include <networkit/community/QuasiThresholdMover/BucketQueue.hpp>
#include <networkit/community/QuasiThresholdMover/DynamicForest.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {
namespace QuasiThresholdMoving {
class EditingRunner {
public:
    EditingRunner(const Graph &G, QuasiThresholdEditingLocalMover::Initialization initialization,
                  count maxIterations, bool sortPaths, bool randomness, count maxPlateauSize,
                  bool useBucketQueue, std::vector<node> order);

    void runLocalMover();

    count getNumberOfEdits() const { return numEdits; };

    count getUsedIterations() const { return usedIterations; };

    count getPlateauSize() const { return actualMaximumPlateau; };

    count getRootEqualBestParents() const { return rootEqualBestParentsCpy; };

    std::map<std::string, std::vector<count>> getRunningInfo() const { return runningInfo; }

    Graph getQuasiThresholdGraph() const;

private:
    struct TraversalData {
        count generation;
        count scoreMax;
        count childCloseness;
        node bestParentBelow;
        /**
         * Logarithm of the number of equally good choices for this node.
         * The logarithm is used as the number of choices can be exponential in the number of nodes.
         * Therefore, the logarithm is guaranteed to be linear in the number of nodes.
         */
        double logEqualBestChoices;
        node numIndifferentChildren;

        TraversalData()
            : generation(none), scoreMax(0), childCloseness(0), bestParentBelow(none),
              logEqualBestChoices(-std::numeric_limits<double>::infinity()), numIndifferentChildren(0){};

        void initialize(count currentGeneration) {
            if (currentGeneration != generation) {
                generation = currentGeneration;
                scoreMax = 0;
                childCloseness = 0;
                bestParentBelow = none;
                logEqualBestChoices = -std::numeric_limits<double>::infinity();
                numIndifferentChildren = 0;
            }
        };

        bool hasChoices() {
            return logEqualBestChoices > -std::numeric_limits<double>::infinity();
        }

        void addEqualChoices(count choices) {
            addLogChoices(std::log(choices));
        }

        void addLogChoices(double logChoices) {
            /* This uses the technique explained on
            https://en.wikipedia.org/w/index.php?title=Log_probability&oldid=954092640#Addition_in_log_space
            to avoid issues when std::exp(logEqualBestChoices) cannot
            be represented as a double anymore. */
            if (logChoices > logEqualBestChoices) {
                logEqualBestChoices =
                    logChoices + std::log1p(std::exp(logEqualBestChoices - logChoices));
            } else {
                logEqualBestChoices += std::log1p(std::exp(logChoices - logEqualBestChoices));
            }
        }

        std::string toString() {
            std::stringstream ss;
            ss << "\n";
            ss << "scoreMax: " << scoreMax << "\n";
            ss << "childCloseness: " << childCloseness << "\n";
            ss << "logEqualBestChoices: " << logEqualBestChoices << "\n";
            ss << "bestParentBelow: " << bestParentBelow << "\n";
            return ss.str();
        };
    };

    const Graph &G;
    count maxIterations;
    count usedIterations;
    bool sortPaths;
    bool randomness;
    std::vector<node> order;
    count maxPlateauSize;

    bool insertRun;
    bool useBucketQueue;

    count numEdits;

    Aux::SignalHandler handler;
    DynamicForest dynamicForest;
    bool hasMoved;
    std::vector<bool> marker;

    count level;
    std::vector<node> neighborQueue;
    std::vector<node> currentLevel;
    std::vector<node> nextLevel;

    BucketQueue bucketQueue;
    std::vector<node> neighbors;
    std::vector<node> touchedNodes;
    std::vector<node> lastVisitedDFSNode;
    std::vector<TraversalData> traversalData;
    std::vector<bool> nodeTouched;

    TraversalData rootData;

    count bestEdits;
    count curEdits;
    node curParent;
    std::vector<node> curChildren;
    std::vector<node> bestChildren;

    std::vector<bool> existing;

    count rootEqualBestParentsCpy;

    count editsBefore;
    count currentPlateau;
    count actualMaximumPlateau;

    count maxDepth;

    std::mt19937_64 &gen;
    std::uniform_real_distribution<double> realDist;
    std::uniform_int_distribution<count> intDist;

    Aux::Timer timer;
    std::vector<std::pair<std::string, Aux::PerfEventCountHardware>> event_counters;
    std::map<std::string, std::vector<count>> runningInfo;

    void localMove(node nodeToMove, count generation);
    void processNode(node u, node nodeToMove, count generation);
    void compareWithQuadratic(node nodeToMove, count generation) const;
    count countNumberOfEdits() const;
    count editsIncidentTo(node u) const;

    bool logRandomBool(double logProbability) {
        assert(logProbability <= 0);
        double x = realDist(gen);
        if (x == 0) return logProbability > - std::numeric_limits<double>::infinity();
        return std::log(x) < logProbability;
    }

    bool randomBool(count options, count optionsToConsider = 1) {
        assert(options > 0);
        assert(options >= optionsToConsider);
        count x = intDist(gen, std::uniform_int_distribution<count>::param_type(0, options - 1));
        return x < optionsToConsider;
    };
};
} // namespace QuasiThresholdMoving

} // namespace NetworKit

#endif // EDITINGRUNNER_H
