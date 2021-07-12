/*
 * LouvainMapEquation.cpp
 *
 * Created on: 2019-01-28
 * Author: Armin Wiebigke
 *         Michael Hamann
 *         Lars Gottesbüren
 */
// networkit-format

#include <algorithm>
#include <cassert>

#include <omp.h>

#include <mutex>
#include <networkit/auxiliary/Parallelism.hpp>
#include <networkit/auxiliary/SignalHandling.hpp>
#include <networkit/coarsening/ParallelPartitionCoarsening.hpp>
#include <networkit/community/LouvainMapEquation.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/structures/Partition.hpp>

namespace NetworKit {

LouvainMapEquation::LouvainMapEquation(const Graph &graph, bool hierarchical, count maxIterations,
                                       ParallelizationType parallelizationType)
    : CloneableCommunityDetectionAlgorithm<LouvainMapEquation>(graph),
      parallel(parallelizationType > ParallelizationType::None),
      parallelizationType(parallelizationType), hierarchical(hierarchical),
      maxIterations(maxIterations) {}

void LouvainMapEquation::run() {
    if (hasRun)
        throw std::runtime_error("Algorithm was already run!");

    Aux::SignalHandler handler;

    result = Partition(G->upperNodeIdBound());
    result.setUpperBound(G->numberOfNodes());
    {
        index i = 0;
        G->forNodes([&](node u) { result[u] = i++; });
    }
    if (parallelizationType == ParallelizationType::Synchronous) {
        nextPartition = Partition(G->upperNodeIdBound());
        {
            index i = 0;
            G->forNodes([&](node u) { nextPartition[u] = i++; });
        }
    }
    if (parallelizationType == ParallelizationType::RelaxMap) {
        locks = Aux::SpinlockArray(G->upperNodeIdBound());
    }
    ets_neighborClusterWeights.resize(parallel ? Aux::getMaxNumberOfThreads() : 1);

    // remove unused nodes from clustering
    if (G->numberOfNodes() != G->upperNodeIdBound()) {
        for (node u = 0; u < G->upperNodeIdBound(); ++u) {
            if (!G->hasNode(u)) {
                result.remove(u);
                if (parallelizationType == ParallelizationType::Synchronous) {
                    nextPartition.remove(u);
                }
            }
        }
    }
    handler.assureRunning();

    calculateInitialClusterCutAndVolume();

#ifndef NDEBUG
    updatePLogPSums();
#endif

    std::vector<AtomicRound> lastActiveRound; // not a member to retain copy constructor

    bool clusteringChanged = false;
    std::vector<node> nodes{G->nodeRange().begin(), G->nodeRange().end()};
    count numberOfNodesMoved = 1;
    for (count iteration = 0; iteration < maxIterations && numberOfNodesMoved > 0; ++iteration) {
        handler.assureRunning();
        std::shuffle(nodes.begin(), nodes.end(), Aux::Random::getURNG());
        if (parallelizationType == ParallelizationType::Synchronous) {
            numberOfNodesMoved = synchronousLocalMoving(nodes, iteration);
        } else if (parallelizationType == ParallelizationType::RelaxMap) {
            numberOfNodesMoved = parallelLocalMoving(nodes, iteration, lastActiveRound);
        } else {
            numberOfNodesMoved = sequentialLocalMoving(nodes, iteration);
        }
        clusteringChanged |= numberOfNodesMoved > 0;

        DEBUG("iteration ", iteration, " num moves ", numberOfNodesMoved);
#ifndef NDEBUG
        DEBUG("Map equation is ", mapEquation());
#endif
    }
    useActiveNodes = false;

    handler.assureRunning();
    if (hierarchical && clusteringChanged) {
        runHierarchical();
    }
    hasRun = true;
}

count LouvainMapEquation::sequentialLocalMoving(std::vector<node> &nodes, count iteration) {
    std::vector<Move> dummyMoves;
    count nodesMoved = 0;
    SparseVector<double> &neighborClusterWeights = ets_neighborClusterWeights[0];
    if (iteration == 0) {
        neighborClusterWeights.resize(result.upperBound(), 0.0);
    }
    for (node u : nodes) {
        if (tryLocalMove(u, neighborClusterWeights, dummyMoves, false)) {
            nodesMoved += 1;
        }
    }
    return nodesMoved;
}

count LouvainMapEquation::parallelLocalMoving(std::vector<node> &nodes, count iteration,
                                              std::vector<AtomicRound> &lastActiveRound) {
    // dummies, since SLM implementation needs more datastructures
    std::vector<Move> dummyMoves;
    count nodesMoved = 0;

    std::mutex activeNodesLock;
    activeNodes.clear();
    auto activateNode = [&](const node u) -> bool {
        auto l = lastActiveRound[u].load(std::memory_order_relaxed);
        if (l != iteration)
            return false;
        return lastActiveRound[u].compare_exchange_strong(l, iteration, std::memory_order_acq_rel);
    };

#pragma omp parallel
    {
        count nm = 0;
        int tid = omp_get_thread_num();
        SparseVector<double> &neighborClusterWeights = ets_neighborClusterWeights[tid];
        if (iteration == 0) {
            neighborClusterWeights.resize(result.upperBound(), 0.0);
        }

        std::vector<node> activeNodesBuffer;

#pragma omp for schedule(guided) nowait
        for (omp_index i = 0; i < static_cast<omp_index>(nodes.size()); ++i) {
            const node u = nodes[i];
            if (tryLocalMove(u, neighborClusterWeights, dummyMoves, false)) {
                nm += 1;
                if (useActiveNodes) {
                    G->forNeighborsOf(u, [&](node v) {
                        if (activateNode(v)) {
                            activeNodesBuffer.push_back(v);
                        }
                    });
                    if (activateNode(u)) {
                        activeNodesBuffer.push_back(u);
                    }
                    // flush local buffer but don't insist on it
                    if (activeNodesBuffer.size() > maximumActiveNodesBufferSize
                        && activeNodesLock.try_lock()) {
                        activeNodes.insert(activeNodes.end(), activeNodesBuffer.begin(),
                                           activeNodesBuffer.end());
                        activeNodesLock.unlock();
                        activeNodesBuffer.clear();
                    }
                }
            }
        }

        // clear out local buffer but insist on the lock.
        if (!activeNodesBuffer.empty()) {
            activeNodesLock.lock();
            activeNodes.insert(activeNodes.end(), activeNodesBuffer.begin(),
                               activeNodesBuffer.end());
            activeNodesLock.unlock();
            activeNodesBuffer.clear();
        }

#pragma omp atomic
        nodesMoved += nm;
    }

    if (useActiveNodes) {
        std::swap(nodes, activeNodes);
    }

    if (nodesMoved > 0 && !useActiveNodes && parallelizationType == ParallelizationType::RelaxMap
        && iteration >= earliestActiveNodesIteration
        && nodesMoved < minimumMoveFraction * G->numberOfNodes()) {
        DEBUG("start using active nodes for the level. threshold = ",
              minimumMoveFraction * G->numberOfNodes(), " num moves = ", nodesMoved);
        useActiveNodes = true;
        if (lastActiveRound.empty()) {
            // assign because std::atomic has no copy constructor
            lastActiveRound = std::vector<AtomicRound>(G->upperNodeIdBound());
        }
    }

    if (useActiveNodes && nodesMoved == 0) {
        // do one full round at the end
        useActiveNodes = false;
        nodes = std::vector<node>(G->nodeRange().begin(), G->nodeRange().end());
    }

    return nodesMoved;
}

count LouvainMapEquation::synchronousLocalMoving(std::vector<node> &nodes, count iteration) {
    // we want at least 10000 nodes per thread in each subround
    count desiredSubroundSize = 10000 * Aux::getCurrentNumberOfThreads();
    // and at least 5 subrounds
    if (nodes.size() >= 5 && nodes.size() / 5 < desiredSubroundSize) {
        desiredSubroundSize = nodes.size() / 5;
    }
    count numSubrounds = idivCeil(nodes.size(), desiredSubroundSize);
    // but at most 64
    if (numSubrounds > 64) {
        numSubrounds = 64;
    }
    const count subroundSize = idivCeil(nodes.size(), numSubrounds);

    count numberOfNodesMoved = 0;

#pragma omp parallel
    {

        int tid = omp_get_thread_num();

        SparseVector<double> &neighborClusterWeights = ets_neighborClusterWeights[tid];
        std::vector<Move> moves;

        // with thread pinning, this ensures the datastructures are allocated on the right socket
        if (iteration == 0) {
            neighborClusterWeights.resize(G->upperNodeIdBound(), 0.0);
        }

        for (size_t i = 0; i < numSubrounds; ++i) {
            size_t first, last;
            std::tie(first, last) = chunkBounds(i, nodes.size(), subroundSize);

            // find moves but don't apply them yet
#pragma omp for schedule(guided)
            for (omp_index j = first; j < last; ++j) {
                const node u = nodes[j];
                tryLocalMove(u, neighborClusterWeights, moves, true);
            }
            // implicit barrier at the end of the loop

            // apply the moves and update cut/volume
            // every thread applies the moves that it calculated
            aggregateAndApplyCutAndVolumeUpdates(moves);
#pragma omp barrier

            // apply moves
            for (Move &m : moves) {
                assert(m.targetCluster == nextPartition[m.movedNode]);
                result[m.movedNode] = m.targetCluster;
            }
#pragma omp atomic
            numberOfNodesMoved += moves.size();
            moves.clear();
#pragma omp barrier

#ifndef NDEBUG
#pragma omp single
            { checkUpdatedCutsAndVolumesAgainstRecomputation(); }
// don't start the next round of moves, before the old ones were validated
#pragma omp barrier
#endif
        }
    }
    return numberOfNodesMoved;
}

void LouvainMapEquation::aggregateAndApplyCutAndVolumeUpdates(std::vector<Move> &moves) {
    double totalCutUpdate = 0.0;
    for (Move &move : moves) {
        const index originCluster = move.originCluster;
        const index targetCluster = move.targetCluster;
        const node u = move.movedNode;

// apply volume updates
#pragma omp atomic
        clusterVolume[originCluster] -= move.volume;
#pragma omp atomic
        clusterVolume[targetCluster] += move.volume;

        double originClusterCutUpdate = move.cutUpdateToOriginCluster,
               targetClusterCutUpdate = move.cutUpdateToTargetCluster;

        // correct the cut for potentially moved neighbors
        // the already calculated cut updates assumed none of them were moved
        G->forEdgesOf(u, [&](node, node v, edgeweight weight) {
            if (u >= v)
                return;

            const index originClusterOfNeighbor = result[v],
                        targetClusterOfNeighbor = nextPartition[v];
            if (originClusterOfNeighbor != targetClusterOfNeighbor) {
                const double w2 = 2 * weight;

                if (originCluster == originClusterOfNeighbor) {
                    originClusterCutUpdate -= w2;
                } else if (targetClusterOfNeighbor == originCluster) {
                    originClusterCutUpdate += w2;
                }

                if (targetClusterOfNeighbor == targetCluster) {
                    targetClusterCutUpdate -= w2;
                } else if (originClusterOfNeighbor == targetCluster) {
                    targetClusterCutUpdate += w2;
                }
            }
        });

// apply cut updates
#pragma omp atomic
        clusterCut[originCluster] += originClusterCutUpdate;
#pragma omp atomic
        clusterCut[targetCluster] += targetClusterCutUpdate;

        totalCutUpdate += originClusterCutUpdate + targetClusterCutUpdate;
    }

#pragma omp atomic
    totalCut += totalCutUpdate;
}

// for every node. store its neighbors that are in the current chunk, and their old cluster IDs, and
// edge weights
bool LouvainMapEquation::tryLocalMove(node u, SparseVector<double> &neighborClusterWeights,
                                      /* SLM specifics */
                                      std::vector<Move> &moves, bool synchronous) {
    // find neighbor clusters
    double vol = 0;
    double loop = 0;
    double weightToCurrent = 0;
    const index currentCluster = result[u];

    G->forEdgesOf(u, [&](node, node v, edgeweight weight) {
        vol += weight;
        if (u != v) {
            const index neighborCluster = result[v];
            if (neighborCluster == currentCluster) {
                weightToCurrent += weight;
            } else {
                if (!neighborClusterWeights.indexIsUsed(neighborCluster))
                    neighborClusterWeights.insert(neighborCluster, 0);
                neighborClusterWeights[neighborCluster] += weight;
            }
        } else {
            loop += weight;
            vol += weight;
        }
    });

    assert(vol == G->weightedDegree(u, true));

    if (neighborClusterWeights.size() > 0) {
        // Calculate best cluster
        index targetCluster = currentCluster;
        double weightToTargetCluster = weightToCurrent;
        double bestChange = fitnessChange(u, vol, loop, currentCluster, currentCluster,
                                          weightToCurrent, weightToCurrent, totalCut);

        neighborClusterWeights.forElements([&](index neighborCluster,
                                               double neighborClusterWeight) {
            const double change = fitnessChange(u, vol, loop, currentCluster, neighborCluster,
                                                neighborClusterWeight, weightToCurrent, totalCut);
            if (change < bestChange
                || (change == bestChange && neighborCluster < targetCluster
                    && targetCluster != currentCluster)) {
                bestChange = change;
                targetCluster = neighborCluster;
                weightToTargetCluster = neighborClusterWeight;
            }
        });

        neighborClusterWeights.reset();

        if (targetCluster != currentCluster) {
            if /* constexpr */ (synchronous) {
                // save move
                const double cutUpdateToCurrentCluster = -vol + 2 * weightToCurrent + 2 * loop;
                const double cutUpdateToTargetCluster = vol - 2 * weightToTargetCluster - 2 * loop;
                moves.emplace_back(u, vol, currentCluster, targetCluster, cutUpdateToCurrentCluster,
                                   cutUpdateToTargetCluster);
                nextPartition[u] = targetCluster;
                return true;
            } else {
                // perform move directly
                return performMove(u, vol, loop, currentCluster, targetCluster,
                                   weightToTargetCluster, weightToCurrent);
            }
        }
    }
    return false;
}

double LouvainMapEquation::fitnessChange(node, double degree, double loopWeight,
                                         node currentCluster, node targetCluster,
                                         double weightToTarget, double weightToCurrent,
                                         double totalCutCurrently) {
    const double cutTarget = clusterCut[targetCluster];
    const double volTarget = clusterVolume[targetCluster];
    const double cutDifferenceCurrent = 2 * weightToCurrent - degree + 2 * loopWeight;

    double values[5]; // allocate values as array for simd loop below

    // Define meaningful aliases for the values
    double &totalCutNew = values[0], &targetClusterCutNew = values[1],
           &targetClusterCutCurrent = values[2], &targetCutPlusVolumeNew = values[3],
           &targetCutPlusVolumeCurrent = values[4];

    if (currentCluster != targetCluster) {
        double cutDifferenceTarget = degree - 2 * weightToTarget - 2 * loopWeight;

        totalCutNew = totalCutCurrently + cutDifferenceCurrent + cutDifferenceTarget;
        targetClusterCutNew = cutTarget + cutDifferenceTarget;
        targetClusterCutCurrent = cutTarget;
        targetCutPlusVolumeNew = cutTarget + cutDifferenceTarget + volTarget + degree;
        targetCutPlusVolumeCurrent = cutTarget + volTarget;
    } else {
        totalCutNew = totalCutCurrently;
        targetClusterCutNew = cutTarget;
        targetClusterCutCurrent = cutTarget + cutDifferenceCurrent;
        targetCutPlusVolumeNew = cutTarget + volTarget;
        targetCutPlusVolumeCurrent = cutTarget + cutDifferenceCurrent + volTarget - degree;
    }

#ifndef _MSC_VER // MSVC doesn't support OpenMP simd, at least with the flags we use
#pragma omp simd
#endif
    for (size_t i = 0; i < 5; ++i) {
        values[i] /= totalVolume;
        if (values[i] > 0) {
            values[i] *= std::log(values[i]);
        } else {
            values[i] = 0;
        }
    }

    return totalCutNew
           + ((targetCutPlusVolumeNew - targetCutPlusVolumeCurrent)
              - (2 * (targetClusterCutNew - targetClusterCutCurrent)));
}

bool LouvainMapEquation::performMove(node u, double degree, double loopWeight, node currentCluster,
                                     node targetCluster, double weightToTarget,
                                     double weightToCurrent) {
    bool moved = true;
    if (parallel) {
        assert(parallelizationType == ParallelizationType::RelaxMap);

        // lock currentCluster and targetCluster
        locks[std::min(currentCluster, targetCluster)].lock();
        locks[std::max(currentCluster, targetCluster)].lock();

        // recompute weightToCurrent and weightToTarget
        weightToCurrent = 0;
        weightToTarget = 0;
        G->forEdgesOf(u, [&](node, node v, edgeweight weight) {
            if (u != v) {
                if (result[v] == currentCluster) {
                    weightToCurrent += weight;
                } else if (result[v] == targetCluster) {
                    weightToTarget += weight;
                }
            }
        });

        const double totalCutCurrently = totalCut;
        const double fitnessCurrent =
            fitnessChange(u, degree, loopWeight, currentCluster, currentCluster, weightToCurrent,
                          weightToCurrent, totalCutCurrently);
        const double fitnessTarget =
            fitnessChange(u, degree, loopWeight, currentCluster, targetCluster, weightToTarget,
                          weightToCurrent, totalCutCurrently);
        if (fitnessTarget >= fitnessCurrent) {
            moved = false;
        }
    }

    if (moved) {
        double cutDifferenceCurrent = 2 * weightToCurrent - degree + 2 * loopWeight;
        double cutDifferenceTarget = degree - 2 * weightToTarget - 2 * loopWeight;
        clusterCut[currentCluster] += cutDifferenceCurrent;
        clusterCut[targetCluster] += cutDifferenceTarget;
        clusterVolume[currentCluster] -= degree;
        clusterVolume[targetCluster] += degree;
        result.moveToSubset(targetCluster, u);

// still protected with atomic because locks are only for clusters
#pragma omp atomic
        totalCut += cutDifferenceCurrent + cutDifferenceTarget;
    }

    if (parallel) {
        // unlock clusters again
        locks[std::max(currentCluster, targetCluster)].unlock();
        locks[std::min(currentCluster, targetCluster)].unlock();
    }

    return moved;
}

void LouvainMapEquation::runHierarchical() {
    assert(result.numberOfSubsets() < result.numberOfElements());
    // free some memory
    clusterVolume.clear();
    clusterVolume.shrink_to_fit();
    clusterCut.clear();
    clusterCut.shrink_to_fit();

    ParallelPartitionCoarsening coarsening(*G, result, parallel && G->numberOfNodes() > 1000000);
    coarsening.run();
    const Graph &metaGraph = coarsening.getCoarseGraph();
    const auto &fineToCoarseMapping = coarsening.getFineToCoarseNodeMapping();

    assert(metaGraph.numberOfNodes() < G->numberOfNodes());

    DEBUG("Run hierarchical with ", metaGraph.numberOfNodes(), " clusters (from ",
          G->numberOfNodes(), " nodes)");

    ParallelizationType para =
        metaGraph.numberOfNodes() > 10000 ? parallelizationType : ParallelizationType::None;
    LouvainMapEquation recursion(metaGraph, true, maxIterations, para);
    recursion.run();
    const Partition &metaPartition = recursion.getPartition();

    G->forNodes([&](node u) { result[u] = metaPartition[fineToCoarseMapping[u]]; });
}

void LouvainMapEquation::calculateInitialClusterCutAndVolume() {
    totalCut = 0.0;
    totalVolume = 0.0;
    clusterCut.assign(G->upperNodeIdBound(), 0.0);
    clusterVolume.assign(G->upperNodeIdBound(), 0.0);

    if (parallel) {
#pragma omp parallel if (G->upperNodeIdBound() > 50000)
        {
            double tCut = 0, tVol = 0;
#pragma omp for schedule(guided)
            for (omp_index u = 0; u < static_cast<omp_index>(G->upperNodeIdBound()); ++u) {
                if (G->hasNode(u)) {
                    G->forEdgesOf(u, [&](node, node v, edgeweight ew) {
                        if (u != v) {
                            clusterCut[u] += ew;
                        } else {
                            ew *= 2;
                        }
                        clusterVolume[u] += ew;
                    });
                }
                tCut += clusterCut[u];
                tVol += clusterVolume[u];
            }

#pragma omp atomic
            totalCut += tCut;
#pragma omp atomic
            totalVolume += tVol;
        }
    } else {
        for (node u = 0; u < G->upperNodeIdBound(); ++u) {
            if (G->hasNode(u)) {
                G->forEdgesOf(u, [&](node, node v, edgeweight ew) {
                    if (u != v) {
                        clusterCut[u] += ew;
                    } else {
                        ew *= 2;
                    }
                    clusterVolume[u] += ew;
                });
            }
            totalCut += clusterCut[u];
            totalVolume += clusterVolume[u];
        }
    }
}

std::string LouvainMapEquation::toString() const {
    return "LouvainMapEquation";
}

#ifndef NDEBUG

double LouvainMapEquation::plogpRel(double w) {
    if (w > 0) {
        double p = w / totalVolume;
        return p * log(p);
    }
    return 0;
}

void LouvainMapEquation::updatePLogPSums() {
    sumPLogPClusterCutPlusVol = 0;
    sumPLogPClusterCut = 0;
    sumPLogPwAlpha = 0;
    for (index i = 0; i < clusterCut.size(); ++i) {
        sumPLogPClusterCutPlusVol += plogpRel(clusterCut[i] + clusterVolume[i]);
        sumPLogPClusterCut += plogpRel(clusterCut[i]);
        sumPLogPwAlpha += plogpRel(clusterVolume[i]);
    }
}

long double LouvainMapEquation::mapEquation() {
    return plogpRel(totalCut) - 2 * sumPLogPClusterCut + sumPLogPClusterCutPlusVol - sumPLogPwAlpha;
}

void LouvainMapEquation::checkUpdatedCutsAndVolumesAgainstRecomputation() {
    std::vector<double> cut(G->upperNodeIdBound(), 0.0), vol(G->upperNodeIdBound(), 0.0);
    double tCut = 0.0, tVol = 0.0;
    G->forNodes([&](node u) {
        double volU = 0.0, cutU = 0.0;
        const index cu = result[u];
        assert(cu == nextPartition[u]);
        G->forEdgesOf(u, [&](node, node v, edgeweight weight) {
            if (cu != result[v]) {
                cutU += weight;
            }

            if (u == v)
                weight *= 2;
            volU += weight;
        });

        vol[result[u]] += volU;
        tVol += volU;

        cut[result[u]] += cutU;
        tCut += cutU;
    });

    G->forNodes([&](node u) {
        const index cu = result[u];
        assert(vol[cu] == clusterVolume[cu]);
        assert(cut[cu] == clusterCut[cu]);
    });
    assert(tVol == totalVolume);
    assert(tCut == totalCut);
}

#endif

} // namespace NetworKit
