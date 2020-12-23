/*
 *
 */
#include <tlx/unused.hpp>
#include <networkit/community/QuasiThresholdEditingLocalMover.hpp>
#include <networkit/community/QuasiThresholdMover/EditingRunner.hpp>
#include <networkit/community/QuasiThresholdMover/QuasiThresholdEditingLinear.hpp>
#include <networkit/generators/TreeReachabilityGraphGenerator.hpp>

namespace NetworKit {

namespace QuasiThresholdMoving {
EditingRunner::EditingRunner(const Graph &G,
                             QuasiThresholdEditingLocalMover::Initialization initialization,
                             count maxIterations, bool sortPaths, bool randomness,
                             count maxPlateauSize, bool useBucketQueue, std::vector<node> order)
    : G(G), maxIterations(maxIterations), usedIterations(0), sortPaths(sortPaths),
      randomness(randomness), maxPlateauSize(maxPlateauSize),
      insertRun(initialization != QuasiThresholdEditingLocalMover::TRIVIAL
                && initialization != QuasiThresholdEditingLocalMover::EDITING),
      useBucketQueue(useBucketQueue), handler(), hasMoved(true),
      marker(G.upperNodeIdBound(), false), lastVisitedDFSNode(G.upperNodeIdBound(), none),
      traversalData(G.upperNodeIdBound()), nodeTouched(G.upperNodeIdBound(), false), rootData(),
      existing(G.upperNodeIdBound(), !insertRun), rootEqualBestParentsCpy(0), currentPlateau(0),
      actualMaximumPlateau(0), gen(Aux::Random::getURNG()), realDist(), intDist(0, 1) {

    runningInfo["time"] = std::vector<count>();
    runningInfo["edits"] = std::vector<count>();
    runningInfo["nodes_moved"] = std::vector<count>();

    if (Aux::PerfEventCountHardware::is_available) {
        event_counters.emplace_back("cache_misses", Aux::PerfEventCountHardware::CACHE_MISSES);
        event_counters.emplace_back("cache_references",
                                    Aux::PerfEventCountHardware::CACHE_REFERENCES);
        event_counters.emplace_back("cycles", Aux::PerfEventCountHardware::CPU_CYCLES);
        event_counters.emplace_back("instructions", Aux::PerfEventCountHardware::INSTRUCTIONS);
    }

    for (const auto &ec : event_counters) {
        runningInfo[ec.first] = std::vector<count>();
    }

    for (auto &ec : event_counters) {
        ec.second.enable();
    }

    timer.start();
    switch (initialization) {
    case QuasiThresholdEditingLocalMover::TRIVIAL: {
        if (G.upperNodeIdBound() == G.numberOfNodes()) {
            dynamicForest = DynamicForest(std::vector<node>(G.upperNodeIdBound(), none));
        } else {
            dynamicForest = DynamicForest(G, std::vector<node>(G.upperNodeIdBound(), none));
        }
        break;
    }
    case QuasiThresholdEditingLocalMover::EDITING: {
        QuasiThresholdEditingLinear editing(G);
        editing.run();
        if (G.upperNodeIdBound() == G.numberOfNodes()) {
            dynamicForest = DynamicForest(editing.getParents());
        } else {
            dynamicForest = DynamicForest(G, editing.getParents());
        }
        break;
    }
    case QuasiThresholdEditingLocalMover::RANDOM_INSERT: {
        this->order.reserve(G.numberOfNodes());
        dynamicForest = DynamicForest(std::vector<node>(G.upperNodeIdBound(), none));
        G.forNodesInRandomOrder([&](node u) { this->order.push_back(u); });
        break;
    }
    case QuasiThresholdEditingLocalMover::ASC_DEGREE_INSERT: {
        dynamicForest = DynamicForest(std::vector<node>(G.upperNodeIdBound(), none));

        std::vector<std::vector<node>> buckets(G.numberOfNodes());
        G.forNodes([&](node u) { buckets[G.degree(u)].push_back(u); });
        this->order.reserve(G.numberOfNodes());
        for (const std::vector<node> &bucket : buckets) {
            for (node u : bucket) {
                this->order.push_back(u);
            }
        }
        break;
    }
    case QuasiThresholdEditingLocalMover::DESC_DEGREE_INSERT: {
        dynamicForest = DynamicForest(std::vector<node>(G.upperNodeIdBound(), none));

        std::vector<std::vector<node>> buckets(G.numberOfNodes());
        G.forNodes([&](node u) { buckets[G.degree(u)].push_back(u); });
        this->order.reserve(G.numberOfNodes());
        for (auto it = buckets.rbegin(); it != buckets.rend(); ++it) {
            for (node u : *it) {
                this->order.push_back(u);
            }
        }
        break;
    }
    case QuasiThresholdEditingLocalMover::USER_DEFINED_INSERT: {
        dynamicForest = DynamicForest(std::vector<node>(G.upperNodeIdBound(), none));
        this->order = std::move(order);
    }
    default:
        break;
    }

    handler.assureRunning();
    if (useBucketQueue) {
        bucketQueue = BucketQueue(G.upperNodeIdBound());
    } else {
        level = 0;
    }

    numEdits = countNumberOfEdits();
    editsBefore = numEdits;

    G.forNodes([&](node u) { lastVisitedDFSNode[u] = u; });

    timer.stop();

    for (auto &ec : event_counters) {
        ec.second.disable();
    }

    for (auto &ec : event_counters) {
        runningInfo[ec.first].push_back(ec.second.readValue());
    }

    runningInfo["time"].push_back(timer.elapsedMicroseconds());
    runningInfo["nodes_moved"].push_back(0);
}

void EditingRunner::runLocalMover() {
    handler.assureRunning();
    if (!insertRun) {
        runningInfo["edits"].push_back(numEdits);
    }
    count generation = 0;
    for (count i = insertRun ? 0 : 1; hasMoved && i <= maxIterations; ++i) {
        if (!hasMoved || (randomness && (currentPlateau >= maxPlateauSize)))
            break;
        handler.assureRunning();
        hasMoved = false;
        numNodesMoved = 0;
        timer.start();

        for (auto &ec : event_counters) {
            ec.second.reset();
            ec.second.enable();
        }

        if (insertRun) {
            for (index j = 0; j < G.numberOfNodes(); j++) {
                node nodeToMove = order[j];
                localMove(nodeToMove, generation++);
                existing[nodeToMove] = 1;
            }
            insertRun = 0;
        } else {
            G.forNodesInRandomOrder([&](node nodeToMove) { localMove(nodeToMove, generation++); });
            INFO("Iteration: ", i, " edits: ", numEdits, " moved nodes: ", numNodesMoved);
        }

        timer.stop();
        for (auto &ec : event_counters) {
            ec.second.disable();
        }

        if (i == 0) {
            runningInfo["nodes_moved"][0] = numNodesMoved;
            runningInfo["time"][0] += timer.elapsedMicroseconds();
            for (auto &ec : event_counters) {
                runningInfo[ec.first][0] += ec.second.readValue();
            }
        } else {
            runningInfo["nodes_moved"].push_back(numNodesMoved);
            runningInfo["time"].push_back(timer.elapsedMicroseconds());
            for (auto &ec : event_counters) {
                runningInfo[ec.first].push_back(ec.second.readValue());
            }
        }
        runningInfo["edits"].push_back(numEdits);
        usedIterations = i;

        assert(numEdits == countNumberOfEdits());
        if (numEdits == editsBefore) {
            currentPlateau++;
        } else {
            if (currentPlateau > actualMaximumPlateau) {
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

void EditingRunner::localMove(node nodeToMove, count generation) {
    assert(numEdits == countNumberOfEdits());
    TRACE("Move node ", nodeToMove);
    handler.assureRunning();
    numNeighbors = 0;
    G.forEdgesOf(nodeToMove, [&](node v) {
        if (!insertRun || existing[v]) {
            ++numNeighbors;
            marker[v] = true;
            neighbors.push_back(v);
        }
    });

    if (!sortPaths) {
        // Do not even attempt to store children when sortPaths is on
        // as due to the sorting they become meaningless.
        curChildren = dynamicForest.children(nodeToMove);
        curParent = dynamicForest.parent(nodeToMove);
    }

    maxDepth = 2 * numNeighbors;
    curEdits = numNeighbors;

    if (!insertRun) {
        // Calculate the old number of edits incident to c to be able to compute the improvement
        // later
        dynamicForest.dfsFrom(
            nodeToMove,
            [&](node c) {
                if (c != nodeToMove) {
                    curEdits += 1 - 2 * marker[c];
                }
            },
            [](node) {});
        dynamicForest.forAncestors(nodeToMove, [&](node p) { curEdits += 1 - 2 * marker[p]; });
    }

    dynamicForest.isolate(nodeToMove);

    if (sortPaths) {
        for (node v : neighbors) {
            dynamicForest.moveUpNeighbor(v, nodeToMove);
        }
    }

    if (useBucketQueue) {
        bucketQueue.fill(neighbors, dynamicForest);
    } else {
        for (node v : neighbors) {
            neighborQueue.emplace_back(v);
        }
        std::stable_sort(neighborQueue.begin(), neighborQueue.end(), [&](node u, node v) {
            return dynamicForest.depth(u) < dynamicForest.depth(v);
        });
    }

    bestChildren.clear();
    rootData.initialize(generation);

    if (useBucketQueue) {
        while (!bucketQueue.empty()) {
            node u = bucketQueue.next();
            processNode(u, nodeToMove, generation);
        }
    } else {
        while (!currentLevel.empty() || !neighborQueue.empty()) {
            if (currentLevel.empty()) {
                level = dynamicForest.depth(neighborQueue.back());
            }
            for (node u : currentLevel) {
                assert(dynamicForest.depth(u) == level);
                processNode(u, nodeToMove, generation);
            }
            while (!neighborQueue.empty() && dynamicForest.depth(neighborQueue.back()) == level) {
                node u = neighborQueue.back();
                neighborQueue.pop_back();
                assert(dynamicForest.depth(u) == level);
                if (nodeTouched[u])
                    continue; // if the node was touched in the previous level, it was in
                              // currentLevel and thus has already been processed
                processNode(u, nodeToMove, generation);
            }
            --level;
            currentLevel.clear();
            currentLevel.swap(nextLevel);
        }
    }

    if (!randomness) {
        if (rootData.childCloseness > rootData.scoreMax) {
            rootData.bestParentBelow = none;
            rootData.scoreMax = rootData.childCloseness;
        }
    } else {
        bool coin = false;
        double ownWeight = rootData.numIndifferentChildren * std::log(2);
        if (rootData.childCloseness > rootData.scoreMax || !rootData.hasChoices()) {
            // INFO("root better");
            rootData.scoreMax = rootData.childCloseness;
            rootData.logEqualBestChoices = ownWeight;
            coin = true;
        } else if (rootData.childCloseness == rootData.scoreMax) {
            ownWeight = rootData.calculateOwnWeightForEqualChoices();
            if (ownWeight > -std::numeric_limits<double>::infinity()) {
                rootData.addLogChoices(ownWeight);
                coin = logRandomBool(ownWeight - rootData.logEqualBestChoices);
            }
            // INFO("root equally good");
        }
        if (coin) {
            rootData.bestParentBelow = none;
        }

        assert(rootData.hasChoices());
    }

    bestEdits = numNeighbors - rootData.scoreMax;

    // If sortPaths and randomness is on, only adopt children when the chosen parent is the
    // lower end of its path.
    const TraversalData &bestParentData =
        (rootData.bestParentBelow == none) ? rootData : traversalData[rootData.bestParentBelow];
    // Check if there are any close children at all, or if there are no close children, if
    // randomness is enabled and we have at least two indifferent children (one alone would not be
    // adopted).
    if (bestParentData.numCloseChildren > 0
        || (randomness && bestParentData.numIndifferentChildren > 1)) {
        std::vector<node> indifferentChildren;
        for (node u : touchedNodes) {
            if (u != nodeToMove && dynamicForest.parent(u) == rootData.bestParentBelow) {
                if (traversalData[u].childCloseness > 0) {
                    bestChildren.push_back(u);
                } else if (randomness && traversalData[u].childCloseness == 0) {
                    indifferentChildren.push_back(u);
                }
            }
        }

        assert(bestChildren.size() == bestParentData.numCloseChildren);

        if (randomness) { // make sure we adopt either 0 or at least two children
            assert(indifferentChildren.size() == bestParentData.numIndifferentChildren);
            // If we shall adopt one child, there must be another indifferent child that we just
            // sample randomly to get at least two
            if (bestChildren.size() == 1) {
                assert(!indifferentChildren.empty());
                index i = Aux::Random::index(indifferentChildren.size());
                bestChildren.push_back(indifferentChildren[i]);
                indifferentChildren[i] = indifferentChildren.back();
                indifferentChildren.pop_back();
            }

            // If there are no best children, sample either 0 or at least two indifferent children
            if (bestChildren.empty() && !indifferentChildren.empty()) {
                assert(indifferentChildren.size() != 1);
                if (indifferentChildren.size() == 2) {
                    // Sample either 0 or two nodes from indifferentChildren
                    if (randomBool(2)) {
                        for (node u : indifferentChildren) {
                            bestChildren.push_back(u);
                        }
                    }
                } else {
                    // sample either 0 or at least two nodes from indifferentChildren
                    std::vector<node> sample;
                    do {
                        sample.clear();
                        for (node u : indifferentChildren) {
                            if (randomBool(2)) {
                                sample.push_back(u);
                            }
                        }
                    } while (sample.size() == 1);

                    for (node u : sample) {
                        bestChildren.push_back(u);
                    }
                }
            } else if (bestChildren.size() > 1) {
                // If there are already two children, just sample randomly from the remaining
                // indifferent children
                for (node u : indifferentChildren) {
                    if (randomBool(2)) {
                        bestChildren.push_back(u);
                    }
                }
            }
        }
    }

#ifndef NDEBUG
    compareWithQuadratic(nodeToMove, generation);
#endif

    // calculate the number of saved edits as comparing the absolute number of edits doesn't make
    // sense
    count savedEdits = curEdits - bestEdits;

    // cleanup for linear move
    for (node u : touchedNodes) {
        lastVisitedDFSNode[u] = u;
        nodeTouched[u] = false;
    }

    assert(!randomness || rootData.hasChoices());

    if (rootData.logEqualBestChoices < std::log(std::numeric_limits<long long>::max())) {
        rootEqualBestParentsCpy = std::llround(std::exp(rootData.logEqualBestChoices));
    } else {
        rootEqualBestParentsCpy = std::numeric_limits<count>::max();
    }

    for (node v : neighbors) {
        marker[v] = false;
    }
    neighbors.clear();
    touchedNodes.clear();

    if (sortPaths || savedEdits > 0 || randomness) {
        dynamicForest.moveToPosition(nodeToMove, rootData.bestParentBelow, bestChildren);
        hasMoved |= (savedEdits > 0 || (randomness && rootEqualBestParentsCpy > 1));
        numNodesMoved += (savedEdits > 0 || (randomness && rootEqualBestParentsCpy > 1));
        numEdits -= savedEdits;
#ifndef NDEBUG
        assert(numEdits == countNumberOfEdits());
#endif
    } else {
        dynamicForest.moveToPosition(nodeToMove, curParent, curChildren);
#ifndef NDEBUG
        assert(numEdits == countNumberOfEdits());
#endif
    }
}

void EditingRunner::processNode(node u, node nodeToMove, count generation) {
    TRACE("Process ", u);
    TRACE("Processing node ", u, " of depth ", dynamicForest.depth(u),
          " (node to move: ", nodeToMove, ")");
    TRACE("Parent: ", dynamicForest.parent(u), ", children: ", dynamicForest.children(u));
    assert(u != nodeToMove);
    tlx::unused(nodeToMove);
    if (useBucketQueue) {
        assert(dynamicForest.depth(u) <= maxDepth);
    }
    if (!nodeTouched[u]) {
        nodeTouched[u] = true;
        touchedNodes.emplace_back(u);
        assert(
            marker[u]); // only marked neighbors may be processed without having been touched before
    }

    traversalData[u].initialize(generation);

    int64_t sumPositiveEdits = traversalData[u].childCloseness;
    assert(traversalData[u].childCloseness >= 0);

    traversalData[u].childCloseness += marker[u];
    traversalData[u].childCloseness -=
        1 - marker[u]; // if (marker[u]) { ++traversalData[u].childCloseness; } else {
                       // --traversalData[u].childCloseness; }

    TRACE("Edit difference before descending: ", traversalData[u].childCloseness);

    assert(!marker[u] || traversalData[u].childCloseness > 0);

    if (traversalData[u].childCloseness >= 0) {
        assert(lastVisitedDFSNode[u] == u);

        node c = dynamicForest.nextDFSNodeOnEnter(u, u);

        while (c != u) {

            if (!nodeTouched[c] || traversalData[c].childCloseness < 0) {

                if (traversalData[u].childCloseness == 0 || dynamicForest.depth(c) > maxDepth) {
                    traversalData[u].childCloseness = -1;
                } else {
                    --traversalData[u].childCloseness;
                }

                // advance to the next starting point for the DFS search.
                c = lastVisitedDFSNode[c];

                if (traversalData[u].childCloseness < 0) {
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

    if (!randomness) {
        if (sumPositiveEdits > traversalData[u].scoreMax || traversalData[u].scoreMax == 0) {
            traversalData[u].scoreMax = sumPositiveEdits;
            traversalData[u].bestParentBelow = u;
        }
    } else {
        bool coin = false;
        double ownWeight = traversalData[u].numIndifferentChildren * std::log(2);
        if (sumPositiveEdits > traversalData[u].scoreMax || !traversalData[u].hasChoices()) {
            // INFO(u, " is better count = 1");
            traversalData[u].scoreMax = sumPositiveEdits;
            traversalData[u].logEqualBestChoices = ownWeight;
            // Either we do not adopt children, or we are at the lower end of a path.
            // Otherwise, there must be a node below u that is at least as good.
            assert(!sortPaths || sumPositiveEdits == 0 || dynamicForest.isLowerEnd(u));
            coin = true;
        } else if (sumPositiveEdits == traversalData[u].scoreMax) {
            ownWeight = traversalData[u].calculateOwnWeightForEqualChoices();
            if (ownWeight > -std::numeric_limits<double>::infinity()) {
                traversalData[u].addLogChoices(ownWeight);
                coin = logRandomBool(ownWeight - traversalData[u].logEqualBestChoices);
            }
            assert(traversalData[u].hasChoices());
            // INFO(u, " equally good count = ", traversalData[u].equalBestParents);
        }
        if (coin) {
            traversalData[u].bestParentBelow = u;
        }
    }

    assert(traversalData[u].scoreMax >= 0);

    traversalData[u].scoreMax += marker[u];

    if (traversalData[u].scoreMax > 0) {
        traversalData[u].scoreMax -= 1 - marker[u];
    }
    TRACE("Maximum gain at ", u, ": ", traversalData[u].scoreMax);
    node p = dynamicForest.parent(u);
    TraversalData &parentData = (p == none) ? rootData : traversalData[p];

    parentData.initialize(generation);

    if ((traversalData[u].scoreMax > 0 || traversalData[u].childCloseness > 0) && p != none) {
        if (useBucketQueue) {
            assert(dynamicForest.depth(p) <= maxDepth);
        }
        if (!nodeTouched[p]) {
            nodeTouched[p] = true;
            touchedNodes.push_back(p);
            if (!useBucketQueue) {
                nextLevel.push_back(p);
            } else if (!marker[p]) { // neighbors already in queue
                bucketQueue.insertParent(p);
            }
        }
    }

    if (traversalData[u].scoreMax > parentData.scoreMax) {
        parentData.logEqualBestChoices = traversalData[u].logEqualBestChoices;
        parentData.scoreMax = traversalData[u].scoreMax;
        parentData.bestParentBelow = traversalData[u].bestParentBelow;
        // INFO(u, " better for ", p);
        // INFO("set count to ", parentData.equalBestParents);
    } else if (randomness && traversalData[u].scoreMax == parentData.scoreMax) {
        parentData.addLogChoices(traversalData[u].logEqualBestChoices);
        if (logRandomBool(traversalData[u].logEqualBestChoices - parentData.logEqualBestChoices)) {
            parentData.bestParentBelow = traversalData[u].bestParentBelow;
        }
        // INFO(u, " equally good for ", p);
        // INFO("increase count by ", traversalData[u].equalBestParents, " to ",
        // parentData.equalBestParents);
    }

    if (traversalData[u].childCloseness >= 0) {
        assert(traversalData[u].childCloseness <= traversalData[u].scoreMax);
        parentData.childCloseness += traversalData[u].childCloseness;
        if (traversalData[u].childCloseness == 0) {
            ++parentData.numIndifferentChildren;
        } else {
            ++parentData.numCloseChildren;
        }
    }

    assert(!dynamicForest.children(u).empty() || traversalData[u].childCloseness == 1);
}

void EditingRunner::compareWithQuadratic(node nodeToMove, count generation) const {
    std::vector<int64_t> missingBelow, missingAbove, existingBelow, existingAbove;
    missingBelow.resize(G.upperNodeIdBound(), 0);
    missingAbove.resize(G.upperNodeIdBound(), 0);
    existingBelow.resize(G.upperNodeIdBound(), 0);
    existingAbove.resize(G.upperNodeIdBound(), 0);
    std::vector<bool> usingDeepNeighbors(G.upperNodeIdBound(), false);
    dynamicForest.forChildrenOf(none, [&](node r) {
        if (existing[r]) {
            dynamicForest.dfsFrom(
                r,
                [&](node u) {
                    if (dynamicForest.depth(u) > maxDepth)
                        usingDeepNeighbors[u] = true;
                    if (u != nodeToMove) {
                        missingBelow[u] = missingAbove[u] = 1 - marker[u];
                        existingBelow[u] = existingAbove[u] = marker[u];
                    }
                    node p = dynamicForest.parent(u);
                    if (p != none) {
                        missingAbove[u] += missingAbove[p];
                        existingAbove[u] += existingAbove[p];
                    }
                },
                [&](node u) {
                    node p = dynamicForest.parent(u);
                    if (p != none) {
                        missingBelow[p] += missingBelow[u];
                        existingBelow[p] += existingBelow[u];
                        if (usingDeepNeighbors[u])
                            usingDeepNeighbors[p] = true;
                    }
                });
        }
    });

    assert(missingBelow[nodeToMove] == 0);
    assert(existingBelow[nodeToMove] == 0);

    if (!sortPaths) {
        bool exactValue = true;
        for (node c : curChildren) {
            missingBelow[nodeToMove] += missingBelow[c];
            existingBelow[nodeToMove] += existingBelow[c];
            if (usingDeepNeighbors[c])
                exactValue = false;
        }

        if (curParent != none) {
            missingAbove[nodeToMove] = missingAbove[curParent];
            existingAbove[nodeToMove] = existingAbove[curParent];
            if (usingDeepNeighbors[curParent])
                exactValue = false;
        }
        if (exactValue) {
            assert(curEdits
                   == numNeighbors - existingAbove[nodeToMove] - existingBelow[nodeToMove]
                          + missingAbove[nodeToMove] + missingBelow[nodeToMove]);
        }
    }

    count minEdits = std::numeric_limits<count>::max();
    std::vector<node> minChildren;
    node minParent = curParent;
    G.forNodes([&](node u) {
        if (u == nodeToMove || usingDeepNeighbors[u] || !existing[u])
            return;
        if (existingBelow[u] >= missingBelow[u]
            || (traversalData[u].generation == generation && traversalData[u].childCloseness > 0)) {
            assert(traversalData[u].childCloseness == existingBelow[u] - missingBelow[u]);
        } else if (nodeTouched[u]) {
            assert(traversalData[u].childCloseness < 0);
        }
    });

    G.forNodes([&](node u) {
        if (dynamicForest.children(u).empty() && marker[u] && !usingDeepNeighbors[u]
            && existing[u]) {
            assert(traversalData[u].childCloseness == 1);
        }
    });

    auto tryEditBelow = [&](node p) {
        if (p == nodeToMove)
            return;

        count edits = numNeighbors;
        if (p != none) {
            edits += missingAbove[p];
            edits -= existingAbove[p];
        }

        std::vector<node> children;
        dynamicForest.forChildrenOf(p, [&](node c) {
            if (c == nodeToMove || usingDeepNeighbors[c] || !existing[c])
                return;
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

    dynamicForest.dfsFrom(
        none, [](node) {}, tryEditBelow);
    tryEditBelow(none);

    assert(minEdits >= bestEdits);

    count childClosenessControl = numNeighbors;
    if (rootData.bestParentBelow != none) {
        childClosenessControl -=
            (existingAbove[rootData.bestParentBelow] - missingAbove[rootData.bestParentBelow]);
    }
    for (node u : bestChildren) {
        childClosenessControl -= (existingBelow[u] - missingBelow[u]);
    }
    if (!sortPaths) {
        TRACE("Current edits: ", curEdits, " (with parent ", curParent, " and current children ",
              curChildren, "), minimum edits: ", minEdits);
    } else {
        TRACE("Current edits: ", curEdits, ", minimum edits: ", minEdits);
    }
    TRACE("Quadratic algorithm wants to have new parent ", minParent, " and new children ",
          minChildren);
    TRACE("Linear algorithm wants to have new parent ", rootData.bestParentBelow,
          " and new children ", bestChildren, " edits: ", childClosenessControl);
    assert(minEdits >= childClosenessControl);

    G.forNodes([&](node u) {
        tlx::unused(u);
        TRACE("Node", marker[u] ? "[x]" : "", " ", u, ", parent ", dynamicForest.parent(u));
    });
}

count EditingRunner::countNumberOfEdits() const {
    // count number of edits that are needed with the initial given forest
    count numExistingEdges = 0;
    count numMissingEdges = 0;
    std::vector<bool> marker(G.upperNodeIdBound());

    dynamicForest.forChildrenOf(none, [&](node r) {
        count depth = 0;
        dynamicForest.dfsFrom(
            r,
            [&](node u) { // on enter
                count upperNeighbors = 0;

                G.forNeighborsOf(u, [&](node v) {
                    if (marker[v])
                        ++upperNeighbors;
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
    while (parent != none) {
        visited[parent] = 1;
        if (!G.hasEdge(parent, u))
            edits++;
        parent = dynamicForest.parent(parent);
    }
    dynamicForest.forChildrenOf(u, [&](node c) {
        dynamicForest.dfsFrom(
            c,
            [&](node w) {
                visited[w] = 1;
                if (!G.hasEdge(w, u))
                    edits++;
            },
            [&](node) {});
    });

    G.forNeighborsOf(u, [&](node w) {
        if (!visited[w])
            edits++;
    });

    return edits;
}
} // namespace QuasiThresholdMoving

} // namespace NetworKit
