#ifndef EDITINGRUNNER_H
#define EDITINGRUNNER_H

#include <networkit/auxiliary/PerfEventCountHardware.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/auxiliary/SignalHandling.hpp>
#include <networkit/auxiliary/Timer.hpp>
#include <networkit/community/QuasiThresholdMover/DynamicForest.hpp>
#include <networkit/community/QuasiThresholdMover/BucketQueue.hpp>
#include <networkit/community/QuasiThresholdEditingLocalMover.hpp>
#include <map>

namespace NetworKit {
  namespace QuasiThresholdMoving {
    class EditingRunner {
    public:
      EditingRunner(const Graph& G, QuasiThresholdEditingLocalMover::Initialization initialization, count maxIterations, bool sortPaths, 
        bool randomness, count maxPlateauSize, bool useBucketQueue, std::vector<node> order);
        
        void runLocalMover();
        
        count getNumberOfEdits() const {
          return numEdits;
        };
        
        count getUsedIterations() const {
          return usedIterations;
        };
        
        count getPlateauSize() const {
          return actualMaximumPlateau;
        };
        
        count getRootEqualBestParents() const {
          return rootEqualBestParentsCpy;
        };
        
        std::map<std::string, std::vector<count>> getRunningInfo() const {
          return runningInfo;
        }
        
        Graph getQuasiThresholdGraph() const;
        
        
      private:
        
        
        struct TraversalData{
          count scoreMax;
          count childCloseness;
          count equalBestParents;
          node bestParentBelow;
          
          TraversalData () : scoreMax(0) , childCloseness(0), equalBestParents(0), bestParentBelow(none){};
          
          void reset(){
            scoreMax = 0;
            childCloseness = 0;
            equalBestParents = 0;
            bestParentBelow = none;
          };
          
          std::string toString(){
            std::stringstream ss;
            ss << "\n";
            ss << "scoreMax: " << scoreMax << "\n";
            ss << "childCloseness: " << childCloseness << "\n";
            ss << "equalBestParents: " << equalBestParents << "\n";
            ss << "bestParentBelow: " << bestParentBelow << "\n";
            return ss.str();
          };
        };
        
        const Graph& G;
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
        
        std::mt19937_64& gen;
        std::uniform_int_distribution<count> dist;
        
        Aux::Timer timer;
        std::vector<std::pair<std::string, Aux::PerfEventCountHardware>> event_counters;
        std::map<std::string, std::vector<count>> runningInfo;
        
        void localMove(node nodeToMove);
        void processNode(node u, node nodeToMove);
        void compareWithQuadratic(node nodeToMove) const;
        count countNumberOfEdits() const;
        count editsIncidentTo(node u) const;
        
        bool randomBool(count options, count optionsToConsider = 1) {
          assert(options > 0);
          assert(options >= optionsToConsider);
          count x = dist(gen, std::uniform_int_distribution<count>::param_type(0, options-1));
          return x < optionsToConsider;
        };
        
        
        
        
      };
  } //namespace QuasiThresholdMoving
  


} // namespace NetworKit

#endif // EDITINGRUNNER_H
