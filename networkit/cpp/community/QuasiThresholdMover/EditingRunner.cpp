/*
 *
 */
#include <networkit/community/QuasiThresholdMover/EditingRunner.hpp>
#include <networkit/generators/TreeReachabilityGraphGenerator.hpp>
#include <networkit/community/QuasiThresholdMover/QuasiThresholdEditingLinear.hpp>
#include <networkit/community/QuasiThresholdEditingLocalMover.hpp>

namespace NetworKit {
  
  namespace QuasiThresholdMoving {
    EditingRunner::EditingRunner(
      const Graph& G, 
      QuasiThresholdEditingLocalMover::Initialization initialization, 
      count maxIterations,
      bool sortPaths,
      bool randomness,
      count maxPlateauSize,
      bool useBucketQueue,
      std::vector<node> order) :
      G(G),
      maxIterations(maxIterations),
      usedIterations(0),
      sortPaths(sortPaths),
      randomness(randomness),
      maxPlateauSize(maxPlateauSize),
      useBucketQueue(useBucketQueue),
      handler(),
      hasMoved(true),
      insertRun(initialization != QuasiThresholdEditingLocalMover::TRIVIAL && initialization != QuasiThresholdEditingLocalMover::EDITING),		
      rootData(),
      rootEqualBestParentsCpy(0),
      currentPlateau(0),
      actualMaximumPlateau(0),
      traversalData(G.upperNodeIdBound()),
      nodeTouched(G.upperNodeIdBound(), false),
      existing(G.upperNodeIdBound(), !insertRun),
      marker(G.upperNodeIdBound(), false),
      lastVisitedDFSNode(G.upperNodeIdBound(), none),
      dist(),
      gen(Aux::Random::getURNG()){
        
        runningInfo["time"] = std::vector<count>();
        runningInfo["edits"] = std::vector<count>();
        
        
        timer.start();
        switch(initialization){
          case QuasiThresholdEditingLocalMover::TRIVIAL:
          {
            if(G.upperNodeIdBound() == G.numberOfNodes()){
              dynamicForest = DynamicForest(std::vector<node>(G.upperNodeIdBound(), none));
            } else {
              dynamicForest = DynamicForest(G, std::vector<node>(G.upperNodeIdBound(), none));
            }
            break;
          }
          case QuasiThresholdEditingLocalMover::EDITING:
          {
            QuasiThresholdEditingLinear editing(G);
            editing.run();
            if(G.upperNodeIdBound() == G.numberOfNodes()){
              dynamicForest = DynamicForest(editing.getParents());
            } else {
              dynamicForest = DynamicForest(G, editing.getParents());
            }
            break;
          }
          case QuasiThresholdEditingLocalMover::RANDOM_INSERT:
          {
            this->order.reserve(G.numberOfNodes());
            dynamicForest = DynamicForest(std::vector<node>(G.upperNodeIdBound(), none));
            G.forNodesInRandomOrder([&](node u) {
              this->order.push_back(u);
            });
            break;
          }
          case QuasiThresholdEditingLocalMover::ASC_DEGREE_INSERT:
          {
            dynamicForest = DynamicForest(std::vector<node>(G.upperNodeIdBound(), none));
            
            std::vector<std::vector<node>> buckets (G.numberOfNodes());
            G.forNodes([&](node u) {
              buckets[G.degree(u)].push_back(u);
            });
            this->order.reserve(G.numberOfNodes());
            for(std::vector<node> bucket : buckets){
              for(node u : bucket){
                this->order.push_back(u);
              }
            }
            break;
          }
          case QuasiThresholdEditingLocalMover::USER_DEFINED_INSERT:
          {	
            dynamicForest = DynamicForest(std::vector<node>(G.upperNodeIdBound(), none));
            this->order = order;
          }
          default:
          break;
        }
        timer.stop();
        runningInfo["time"].push_back(timer.elapsedMilliseconds());
        
        
        handler.assureRunning();
        if(useBucketQueue){
          bucketQueue = BucketQueue(G.upperNodeIdBound());
        } else {
          level = 0;
        }
        
        
        numEdits = countNumberOfEdits();
        editsBefore = numEdits;
        
        
        G.forNodes([&](node u) {
          lastVisitedDFSNode[u] = u;
        });
        

        
      }
      
      void EditingRunner::runLocalMover(){
        handler.assureRunning();
        if(!insertRun) {
          runningInfo["edits"].push_back(numEdits);
        }
        for (count i = insertRun ? 0 : 1; hasMoved && i <= maxIterations; ++i) {
          if(!hasMoved || (randomness && (currentPlateau >= maxPlateauSize))) break;
          handler.assureRunning();
          hasMoved = false;
          timer.start();
          if(insertRun){
            for(index j = 0; j < G.numberOfNodes(); j++){
              node nodeToMove = order[j];
              localMove(nodeToMove);
              existing[nodeToMove] = 1;
            }
            insertRun = 0;
          } else {
            G.forNodesInRandomOrder([&](node nodeToMove) {
              localMove(nodeToMove);
            });
          }
          

          timer.stop();
          if(i == 0){
            runningInfo["time"][0] += timer.elapsedMilliseconds();
          } else {
            runningInfo["time"].push_back(timer.elapsedMilliseconds());
          }
          runningInfo["edits"].push_back(numEdits);
          usedIterations = i;
          
          assert(numEdits == countNumberOfEdits());
          if(numEdits == editsBefore){
            currentPlateau++;
          }	else {
            if(currentPlateau > actualMaximumPlateau){
              actualMaximumPlateau = currentPlateau;
            }
            currentPlateau = 0;
          }
          editsBefore = numEdits;
        }
      }
      
      Graph EditingRunner::getQuasiThresholdGraph() const {
        Graph forest = dynamicForest.toGraph();
        TreeReachabilityGraphGenerator gen(forest);
        gen.run();
        return gen.getGraph();
      }
      
      void EditingRunner::localMove(node nodeToMove){
        assert(numEdits == countNumberOfEdits());
        TRACE("Move node ", nodeToMove);
        handler.assureRunning();
        G.forEdgesOf(nodeToMove, [&](node v) {
          if(existing[v]){
            marker[v] = true;
            neighbors.push_back(v);
            if(sortPaths) dynamicForest.moveUpNeighbor(v, nodeToMove);
          }
        });
        if(sortPaths) dynamicForest.moveUpNeighbor(nodeToMove, nodeToMove);
        curChildren = dynamicForest.children(nodeToMove);
        curParent = dynamicForest.parent(nodeToMove);
        if(insertRun){
          maxDepth = 2 * neighbors.size();
          curEdits = neighbors.size();
        } else {
          maxDepth = 2 * neighbors.size();
          curEdits = G.degree(nodeToMove);
          dynamicForest.dfsFrom(nodeToMove, [&](node c) {
            if (c != nodeToMove) {
              curEdits += 1 - 2 * marker[c];
            }
          }, [](node){});
          for (node p = dynamicForest.parent(nodeToMove); p != none; p = dynamicForest.parent(p)) {
            curEdits += 1 - 2 * marker[p];
          }
        }
        dynamicForest.isolate(nodeToMove);
        if(useBucketQueue){
          bucketQueue.fill(neighbors, dynamicForest);
        } else {
          for(node v : neighbors){
            neighborQueue.emplace_back(v);
          }
          std::stable_sort(neighborQueue.begin(), neighborQueue.end(), [&](node u, node v) {return dynamicForest.depth(u) < dynamicForest.depth(v);});
          count level = 0;
          if (!neighborQueue.empty()) {
            level = dynamicForest.depth(neighborQueue.back());
          }
        }
        
        bestChildren.clear();
        rootData.reset();
        //all neighbors to deep
        if(useBucketQueue && bucketQueue.empty() || !useBucketQueue && neighborQueue.empty()) {
          rootData.bestParentBelow = none;
          bestEdits = neighbors.size();
          TRACE("Isolate");
        } else {
          if(useBucketQueue){
            while(!bucketQueue.empty()){
              node u = bucketQueue.next();
              processNode(u, nodeToMove);
            }
          } else {
            while (!currentLevel.empty() || !neighborQueue.empty()) {
              if (currentLevel.empty()) {
                level = dynamicForest.depth(neighborQueue.back());
              }
              for (node u : currentLevel) {
                assert(dynamicForest.depth(u) == level);
                processNode(u, nodeToMove);
              }
              while (!neighborQueue.empty() && dynamicForest.depth(neighborQueue.back()) == level) {
                node u = neighborQueue.back();
                neighborQueue.pop_back();
                assert(dynamicForest.depth(u) == level);
                if (nodeTouched[u]) continue; // if the node was touched in the previous level, it was in currentLevel and thus has already been processed
                processNode(u, nodeToMove);
              }
              --level;
              currentLevel.clear();
              currentLevel.swap(nextLevel);
            }
          }
          
          if(!randomness) {
            if (rootData.childCloseness > rootData.scoreMax) {
              rootData.bestParentBelow = none;
              rootData.scoreMax = rootData.childCloseness;
            } 
          } else {
            bool coin = false;
            if (rootData.childCloseness > rootData.scoreMax || rootData.equalBestParents == 0) {
              //INFO("root better");
              rootData.scoreMax = rootData.childCloseness;
              rootData.equalBestParents = 1;
              coin = true;
            } else if (rootData.childCloseness == rootData.scoreMax) {
              //INFO("root equally good");
              rootData.equalBestParents += 1;
              coin = randomBool(rootData.equalBestParents);
            }
            if (coin) {
              rootData.bestParentBelow = none;
            }
          }
          bestEdits = neighbors.size() - rootData.scoreMax;
          
          
          for (node u : touchedNodes) {
            if (u != nodeToMove && dynamicForest.parent(u) == rootData.bestParentBelow && traversalData[u].childCloseness != none) {
              if(traversalData[u].childCloseness > 0 || (randomness &&  randomBool(2))){
                bestChildren.push_back(u);
              }
            }
          }
        }
        
        #ifndef NDEBUG
        compareWithQuadratic(nodeToMove);
        #endif
        
        // calculate the number of saved edits as comparing the absolute number of edits doesn't make sense
        count savedEdits = curEdits - bestEdits;
        
        // cleanup for linear move
        for (node u : touchedNodes) {
          lastVisitedDFSNode[u] = u;
          nodeTouched[u] = false;
          traversalData[u].reset();
        }
        
        
        neighbors.clear();
        touchedNodes.clear();
        rootEqualBestParentsCpy = rootData.equalBestParents;
        
        
        G.forEdgesOf(nodeToMove, [&](node v) {
          marker[v] = false;
        });
        
        
        if (savedEdits > 0 || (savedEdits == 0 && randomness)) {
          dynamicForest.moveToPosition(nodeToMove, rootData.bestParentBelow, bestChildren);
          hasMoved = true;
          numEdits -= savedEdits;
          #ifndef NDEBUG
          assert(numEdits == countNumberOfEdits());
          #endif
        } else  {
          dynamicForest.moveToPosition(nodeToMove, curParent, curChildren);
          #ifndef NDEBUG
          assert(numEdits == countNumberOfEdits());
          #endif
        }
      }
      
      
      void EditingRunner::processNode(node u, node nodeToMove){
        TRACE("Process ", u);
        TRACE("Processing node ", u, " of depth ", dynamicForest.depth(u), " (node to move: ", nodeToMove, ")");
        TRACE("Parent: ", dynamicForest.parent(u), ", children: ", dynamicForest.children(u));
        assert(u != nodeToMove);
        if(useBucketQueue){
          assert(dynamicForest.depth(u) <= maxDepth);
        }
        if (!nodeTouched[u]) {
          nodeTouched[u] = true;
          touchedNodes.emplace_back(u);
          assert(marker[u]); // only marked neighbors may be processed without having been touched before
        }
        
        count sumPositiveEdits = traversalData[u].childCloseness;
        assert(traversalData[u].childCloseness != none);
        
        traversalData[u].childCloseness += marker[u];
        traversalData[u].childCloseness -= 1 - marker[u]; // if (marker[u]) { ++traversalData[u].childCloseness; } else { --traversalData[u].childCloseness; }
        
        TRACE("Edit difference before descending: ", traversalData[u].childCloseness);
        
        assert(!marker[u] || traversalData[u].childCloseness > 0);
        
        if (traversalData[u].childCloseness != none) {
          assert(lastVisitedDFSNode[u] == u);
          
          node c = dynamicForest.nextDFSNodeOnEnter(u, u);
          
          while (c != u) {
            
            if (!nodeTouched[c] || traversalData[c].childCloseness == none) {
              
              if(dynamicForest.depth(c) > maxDepth){
                traversalData[u].childCloseness = none;
              } else {
                --traversalData[u].childCloseness;
              }
              
              // advance to the next starting point for the DFS search.
              c = lastVisitedDFSNode[c];
              
              if (traversalData[u].childCloseness == none) {
                lastVisitedDFSNode[u] = c;
                break;
              }
              
              c = dynamicForest.nextDFSNodeOnEnter(c, u);
            } else {
              node p = dynamicForest.parent(c);
              c = dynamicForest.nextChild(c, p);
              
              while (c == p && c != u) {
                p = dynamicForest.parent(p);
                c = dynamicForest.nextChild(c, p);
              }
            }
          }
        }
        
        TRACE("Edit difference after descending: ", traversalData[u].childCloseness);
        
        if(!randomness){
          if (sumPositiveEdits > traversalData[u].scoreMax || traversalData[u].scoreMax == 0) {
            traversalData[u].scoreMax = sumPositiveEdits;
            traversalData[u].bestParentBelow = u;
          }
        } else {
          bool coin  = false;
          if (sumPositiveEdits > traversalData[u].scoreMax || traversalData[u].equalBestParents == 0) {
            //INFO(u, " is better count = 1");
            traversalData[u].scoreMax = sumPositiveEdits;
            traversalData[u].equalBestParents = 1;
            coin = true;
          } else if (sumPositiveEdits == traversalData[u].scoreMax) {
            traversalData[u].equalBestParents += 1;
            coin = randomBool(traversalData[u].equalBestParents);
            //INFO(u, " equally good count = ", traversalData[u].equalBestParents);
          }
          if (coin) {
            traversalData[u].bestParentBelow = u;
          }
        }
        
        

        
        assert(traversalData[u].scoreMax != none);
        
        traversalData[u].scoreMax += marker[u];
        
        if (traversalData[u].scoreMax > 0) {
          traversalData[u].scoreMax -= 1 - marker[u];
        }
        TRACE("Maximum gain at ", u, ": ", traversalData[u].scoreMax);
        node p = dynamicForest.parent(u);
        TraversalData parentData;
        if(p == none){
          parentData = rootData;
        } else {
          parentData = traversalData[p];
        }
        if ((traversalData[u].scoreMax > 0 || (traversalData[u].childCloseness != none && traversalData[u].childCloseness > 0)) && p != none){
          if(useBucketQueue){
            assert(dynamicForest.depth(p) <= maxDepth);
          }
          if (!nodeTouched[p]) {
            nodeTouched[p] = true;
            touchedNodes.push_back(p);
            if(!useBucketQueue){
              nextLevel.push_back(p);
            }else if(!marker[p]){  //neighbors already in queue
              bucketQueue.insertParent(p);
            } 
          }
        }
        bool coin = 0;
        if (randomness) {
          if (traversalData[u].scoreMax > parentData.scoreMax){
            parentData.equalBestParents = traversalData[u].equalBestParents;
            //INFO(u, " better for ", p);
            //INFO("set count to ", parentData.equalBestParents);
          }
          if(traversalData[u].scoreMax == parentData.scoreMax){
            if(parentData.scoreMax > 0 || p == none || marker[p]){
              parentData.equalBestParents+=traversalData[u].equalBestParents;
              coin = randomBool(parentData.equalBestParents, traversalData[u].equalBestParents);
              //INFO(u, " equally good for ", p);
              //INFO("increase count by ", traversalData[u].equalBestParents, " to ", parentData.equalBestParents);
            }
          }
        }
        if (traversalData[u].scoreMax > parentData.scoreMax || coin) {
          parentData.scoreMax = traversalData[u].scoreMax;
          parentData.bestParentBelow = traversalData[u].bestParentBelow;
        }
        if (traversalData[u].childCloseness != none) {
          assert(traversalData[u].childCloseness <= traversalData[u].scoreMax);
          parentData.childCloseness += traversalData[u].childCloseness;
        }
        if(p == none){
          rootData = parentData;
        } else {
          traversalData[p] = parentData;
        }
        
        assert(!dynamicForest.children(u).empty() || traversalData[u].childCloseness == 1);
      }
      
      void EditingRunner::compareWithQuadratic(node nodeToMove) const {
        std::vector<count> missingBelow, missingAbove, existingBelow, existingAbove;
        missingBelow.resize(G.upperNodeIdBound(), 0);
        missingAbove.resize(G.upperNodeIdBound(), 0);
        existingBelow.resize(G.upperNodeIdBound(), 0);
        existingAbove.resize(G.upperNodeIdBound(), 0);
        std::vector<bool> usingDeepNeighbors(G.upperNodeIdBound(), false);
        dynamicForest.forChildrenOf(none, [&](node r) {
          if(existing[r]){
            dynamicForest.dfsFrom(r, [&](node u) {
              if(dynamicForest.depth(u) > maxDepth) usingDeepNeighbors[u] = true;
              if (u != nodeToMove) {
                missingBelow[u] = missingAbove[u] = 1 - marker[u];
                existingBelow[u] = existingAbove[u] = marker[u];
              }
              node p = dynamicForest.parent(u);
              if (p != none) {
                missingAbove[u] += missingAbove[p];
                existingAbove[u] += existingAbove[p];
              }
            },	[&](node u) {
              node p = dynamicForest.parent(u);
              if (p != none) {
                missingBelow[p] += missingBelow[u];
                existingBelow[p] += existingBelow[u];
                if(usingDeepNeighbors[u]) usingDeepNeighbors[p] = true;
              }
            });
          }
        });
        
        assert(missingBelow[nodeToMove] == 0);
        assert(existingBelow[nodeToMove] == 0);
        
        bool exactValue = true;
        for (node c : curChildren) {
          missingBelow[nodeToMove] += missingBelow[c];
          existingBelow[nodeToMove] += existingBelow[c];
          if(usingDeepNeighbors[c]) exactValue = false;
        }
        
        if (curParent != none) {
          missingAbove[nodeToMove] = missingAbove[curParent];
          existingAbove[nodeToMove] = existingAbove[curParent];
          if(usingDeepNeighbors[curParent]) exactValue = false;
        }
        if(exactValue){
          assert(curEdits == neighbors.size() - existingAbove[nodeToMove] - existingBelow[nodeToMove] + missingAbove[nodeToMove] + missingBelow[nodeToMove]);
        }
        
        count minEdits = curEdits;
        std::vector<node> minChildren;
        node minParent = curParent;
        G.forNodes([&](node u) {
          if (u == nodeToMove || usingDeepNeighbors[u] || !existing[u]) return;
          if (existingBelow[u] >= missingBelow[u] || (traversalData[u].childCloseness > 0 && traversalData[u].childCloseness != none)) {
            assert(traversalData[u].childCloseness == existingBelow[u] - missingBelow[u]);
          } else if (nodeTouched[u]) {
            assert(traversalData[u].childCloseness == none);
          }
        });
        
        G.forNodes([&](node u) {
          if (dynamicForest.children(u).empty() && marker[u] && !usingDeepNeighbors[u] &&  existing[u]) {
            assert(traversalData[u].childCloseness == 1);
          }
        });
        
        auto tryEditBelow = [&](node p) {
          if (p == nodeToMove) return;
          
          count edits = neighbors.size();
          if (p != none) {
            edits += missingAbove[p];
            edits -= existingAbove[p];
          }
          
          std::vector<node> children;
          dynamicForest.forChildrenOf(p, [&](node c) {
            if (c == nodeToMove || usingDeepNeighbors[c] || !existing[c]) return;
            if (existingBelow[c] > missingBelow[c]) { // TODO try >= (more children...)
              if (dynamicForest.children(c).empty() && marker[c]) {
                assert(traversalData[c].childCloseness == 1);
              }
              assert(traversalData[c].childCloseness == existingBelow[c] - missingBelow[c]);
              
              children.emplace_back(c);
              edits -= existingBelow[c] - missingBelow[c];
            }
          });
          
          if (edits < minEdits) {
            minEdits = edits;
            minChildren = std::move(children);
            minParent = p;
          }
        };
        
        dynamicForest.dfsFrom(none, [](node){}, tryEditBelow);
        tryEditBelow(none);
        
        assert(minEdits >= bestEdits);
        
        count childClosenessControl = neighbors.size();
        if (rootData.bestParentBelow != none) {
          childClosenessControl -= (existingAbove[rootData.bestParentBelow] - missingAbove[rootData.bestParentBelow]);
        }
        for (node u : bestChildren) {
          childClosenessControl -= (existingBelow[u] - missingBelow[u]);
        }
        TRACE("Current edits: ", curEdits, " (with parent ", curParent, " and current children ", curChildren, "), minimum edits: ", minEdits);
        TRACE("Quadratic algorithm wants to have new parent ", minParent, " and new children ", minChildren);
        TRACE("Linear algorithm wants to have new parent ", rootData.bestParentBelow, " and new children ", bestChildren, " edits: ", childClosenessControl);
        assert(minEdits >= childClosenessControl);
      }
      
      count EditingRunner::countNumberOfEdits() const {
        // count number of edits that are needed with the initial given forest
        count numExistingEdges = 0;
        count numMissingEdges = 0;
        std::vector<bool> marker(G.upperNodeIdBound());
        
        dynamicForest.forChildrenOf(none, [&](node r) {
          count depth = 0;
          dynamicForest.dfsFrom(r,
            [&](node u) { // on enter
              count upperNeighbors = 0;
              
              G.forNeighborsOf(u, [&](node v) {
                if (marker[v]) ++upperNeighbors;
              });
              
              numExistingEdges += upperNeighbors;
              numMissingEdges += depth - upperNeighbors;
              marker[u] = true;
              depth += 1;
            },
            [&](node u) { // on exit
              marker[u] = false;
              depth -= 1;
            });
          });
          
          return numMissingEdges + (G.numberOfEdges() - numExistingEdges);
        }
        
        count EditingRunner::editsIncidentTo(node u) const {
          std::vector<bool> visited(G.upperNodeIdBound());
          count edits = 0;
          node parent = dynamicForest.parent(u);
          while(parent != none){
            visited[parent] = 1;
            if(!G.hasEdge(parent, u))	edits++;
            parent = dynamicForest.parent(parent);
          }
          dynamicForest.forChildrenOf(u, [&](node c) {
            dynamicForest.dfsFrom(c,
              [&](node w) {
                visited[w] = 1; 
                if(!G.hasEdge(w, u)) edits++;
              },
              [&](node w) {});
            });
            
            G.forNeighborsOf(u, [&](node w) {
              if (!visited[w]) edits++;
            });
            
            return edits;
          }
  } //namespace QuasiThresholdMoving
  

} // namespace NetworKit


