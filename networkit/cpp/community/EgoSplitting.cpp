/*
 * EgoSplitting.cpp
 *
 * Created: 2018-12-11
 * Author: Armin Wiebigke
 */

#include <cmath>
#include <limits>
#include <stdexcept>
#include <algorithm>

#include <tlx/container.hpp>
#include <networkit/auxiliary/IncrementalUniformRandomSelector.hpp>
#include <networkit/auxiliary/Parallel.hpp>
#include <networkit/auxiliary/ParseString.hpp>
#include <networkit/auxiliary/SignalHandling.hpp>
#include <networkit/auxiliary/VectorComparator.hpp>
#include <networkit/community/ConductanceCommunityCleanup.hpp>
#include <networkit/community/EgoSplitting.hpp>
#include <networkit/community/PLM.hpp>
#include <networkit/community/LouvainMapEquation.hpp>
#include <networkit/community/cleanup/SignificanceCommunityCleanUp.hpp>
#include <networkit/components/ConnectedComponents.hpp>
#include <networkit/coarsening/ParallelPartitionCoarsening.hpp>
#include <networkit/graph/RandomMaximumSpanningForest.hpp>
#include <networkit/structures/Partition.hpp>


namespace NetworKit {

EgoSplitting::EgoSplitting(const Graph &G, bool parallelEgoNetEvaluation)
        : EgoSplitting(G, parallelEgoNetEvaluation,
                       PLM(G, true, 1.0, "none randomized"),
                       LouvainMapEquation(G, true, 64, "relaxmap"))
                       {}

EgoSplitting::EgoSplitting(const Graph &G, bool parallelEgoNetEvaluation,
                           const CommunityDetectionAlgorithm &clusterAlgo)
        : EgoSplitting(G, parallelEgoNetEvaluation, clusterAlgo, clusterAlgo) {
}

EgoSplitting::EgoSplitting(const Graph &G, bool parallelEgoNetEvaluation,
                           const CommunityDetectionAlgorithm &localClusterAlgo,
                           const CommunityDetectionAlgorithm &globalClusterAlgo)
        : Algorithm(), G(&G),
          parallelEgoNetEvaluation(parallelEgoNetEvaluation),
          localClusteringAlgo(localClusterAlgo.clone()),
          globalClusteringAlgo(globalClusterAlgo.clone()),
          stochasticDistribution(0) {
    setDefaultParameters();
}

void EgoSplitting::setDefaultParameters() {
    // Parameters for ego-net extension
    parameters["Extend EgoNet"] = "Yes";
    parameters["Maximum Extend Factor"] = "5";
    parameters["Maximum Extend Exponent"] = "0.5";
    parameters["Minimum Neighbors"] = "3";
    parameters["Extend Only With Conductance Increase"] = "Yes";

    // Connect Personas
    parameters["Connect Personas Internally"] = "Yes";

    // Cleanup
    parameters["Cleanup"] = "Yes";
    parameters["Cleanup Merge"] = "Yes";
    parameters["Cleanup Min Overlap"] = "0.5";
    parameters["Cleanup Conductance"] = "Yes";

    // Debug/Analysis
    parameters["Store EgoNet Graphs"] = "No";
    parameters["Max Stored EgoNets"] = "2000";
}

void EgoSplitting::run() {
    if (hasRun)
        throw std::runtime_error("Algorithm has already been run!");
    if (G->numberOfSelfLoops() > 0)
        throw std::runtime_error("No self-loops allowed!");
    assert(timingsEmpty());
    Aux::SignalHandler handler;
    Aux::Timer timer;
    timer.start();

    INFO("Initializing data structures");
    init();
    addTime(timer, "0  Data Initialization");
    handler.assureRunning();

    INFO("create EgoNets");
    createEgoNets();
    addTime(timer, "1  Create EgoNets");
    handler.assureRunning();

    INFO("split into Personas");
    splitIntoPersonas();
    addTime(timer, "2  Split Personas");
    handler.assureRunning();

    INFO("Connect Personas Internally");
    connectPersonas();
    addTime(timer, "3  Connect Personas");
    handler.assureRunning();

    INFO("create Persona Clustering");
    createPersonaClustering();
    addTime(timer, "4  Persona Clustering");
    handler.assureRunning();

    INFO("create Communities");
    std::vector<std::vector<node>> communities = getCommunitiesFromPersonaClustering();
    addTime(timer, "5  Create Communities");

    personaGraph = Graph(); // deallocate persona graph
    addTime(timer, "6  Dealloc Persona Graph");

    INFO("clean up Communities");
    cleanUpCommunities(communities);
    addTime(timer, "7  Cleanup Communities");

    INFO("create Cover");
    createCover(communities);
    addTime(timer, "8  Create Cover");

    hasRun = true;
}

void EgoSplitting::init() {
    personaEdges.resize(G->upperNodeIdBound());
    egoNetPartitionsOffset.resize(G->upperNodeIdBound() + 1);
    egoNetPartitions.resize(2 * G->numberOfEdges(), std::make_pair(none, none));
    egoNetPartitionCounts.resize(G->upperNodeIdBound(), 0);
    personaOffsets.resize(G->upperNodeIdBound() + 1, 0);
    directedG = LowToHighDirectedGraph(*G);

    if (parameters.at("Cleanup") == "Yes" || parameters.at("Extend EgoNet") == "Yes") {
        stochasticDistribution.increaseMaxValueTo(2 * G->numberOfEdges() + G->numberOfNodes());
    }

    if (parameters.at("Store EgoNet Graphs") == "Yes") {
        egoNetExtendedPartitions.resize(G->upperNodeIdBound());
    }
}

void EgoSplitting::createEgoNets() {
    {
        index sum = 0;
        for (node u = 0; u < egoNetPartitionsOffset.size(); ++u) {
            egoNetPartitionsOffset[u] = sum;
            if (G->hasNode(u)) {
                sum += G->degree(u);
            }
        }
    }
    const bool extendEgoNet = parameters.at("Extend EgoNet") == "Yes";
    const double extendFactor = Aux::stringToDouble(parameters.at("Maximum Extend Factor"));
    const double extendExponent = Aux::stringToDouble(parameters.at("Maximum Extend Exponent"));
    const bool storeEgoNetGraphs = parameters.at("Store EgoNet Graphs") == "Yes";
    const bool connectPersonas = parameters.at("Connect Personas Internally") == "Yes";
    const count minNeighbors = std::stoi(parameters.at("Minimum Neighbors"));
    const bool conductanceExtend = parameters.at("Extend Only With Conductance Increase") == "Yes";
#pragma omp parallel if (parallelEgoNetEvaluation)
    {
        Aux::SignalHandler signalHandler;
        NodeMapping egoMapping(*G); // Assign local IDs to the neighbors
        SparseVector<double> weightToOriginalEgoNet(G->upperNodeIdBound());
        std::vector<std::pair<node, double>> candidatesAndScores;

        auto upperNodeBound = static_cast<omp_index>(G->upperNodeIdBound());
#pragma omp for schedule(dynamic, 10)
        for (omp_index egoNode = 0; egoNode < upperNodeBound; ++egoNode) {
            if (!signalHandler.isRunning()) continue;
            if (!G->hasNode(egoNode) || G->degree(egoNode) < 2) continue;

            egoMapping.reset();

            DEBUG("Create EgoNet for Node ", egoNode, "/", G->upperNodeIdBound());
            // Find neighbors == nodes of the ego-net
            G->forEdgesOf(egoNode, [&](node, node v) {
                if (G->degree(v) > 1) {
                    egoMapping.addNode(v);
                }
            });

            const count originalEgoNetSize = egoMapping.nodeCount();
            if (originalEgoNetSize == 0)
                continue; // egoNode has only degree-1-neighbors

            if (extendEgoNet) {
                // Calculate number of neighbors in egoNet for all neighbors of the egonet
                weightToOriginalEgoNet.reset();
                double egonetVolume = 0;
                double egonetCut = 0;
                for (node u : egoMapping.globalNodes()) {
                    G->forEdgesOf(u, [&](node, node neighbor, edgeweight weight) {
                        egonetVolume += weight;
                        if (neighbor != egoNode && !egoMapping.isMapped(neighbor)) {
                            egonetCut += weight;
                            if (!weightToOriginalEgoNet.indexIsUsed(neighbor)) {
                                weightToOriginalEgoNet.insert(neighbor, weight);
                            } else {
                                weightToOriginalEgoNet[neighbor] += weight;
                            }
                        }
                    });
                }

                // Get scores of extension candidates
                candidatesAndScores.clear();
                weightToOriginalEgoNet.forElements([&](node candidate, double neighborCount) {
                    if (neighborCount >= minNeighbors) { // FIXME: this absolute limit doesn't work for weighted
                        double score = neighborCount * neighborCount / G->degree(candidate); // FIXME: replace degree by (stored) weighted degree for weighted
                        candidatesAndScores.emplace_back(candidate, score);
                    }
                });

                // Select best candidates
                count extendedNodesLimit =
                    std::ceil(extendFactor * std::pow(originalEgoNetSize, extendExponent));
                if (candidatesAndScores.size() > extendedNodesLimit) {
                    std::nth_element(
                        candidatesAndScores.begin(),
                        candidatesAndScores.begin() + extendedNodesLimit, candidatesAndScores.end(),
                        [](const std::pair<node, double> &a, const std::pair<node, double> &b) {
                            return a.second > b.second;
                        });
                    candidatesAndScores.resize(extendedNodesLimit);
                }

                if (conductanceExtend) {
                    // Only add candidates that improve the conductance
                    std::sort(candidatesAndScores.begin(), candidatesAndScores.end(),
                              [](const std::pair<node, double> &a,
                                 const std::pair<node, double> &b) { return a.second > b.second; });

                    for (const std::pair<node, double> &cs : candidatesAndScores) {
                        node candidate = cs.first;
                        double candidateDegree = G->degree(candidate);
                        double edgesToEgoNet = weightToOriginalEgoNet[candidate];

                        double newVolume = egonetVolume + candidateDegree;
                        double newCut = egonetCut - 2 * edgesToEgoNet + candidateDegree;

                        if (egonetCut / egonetVolume >= newCut / newVolume) {
                            egoMapping.addNode(candidate);
                            egonetCut = newCut;
                            egonetVolume = newVolume;
                        }
                    }
                } else {
                    for (const std::pair<node, double> &cs : candidatesAndScores) {
                        egoMapping.addNode(cs.first);
                    }
                }
            }

            Graph egoGraph(egoMapping.nodeCount(), G->isWeighted());

            // Find all triangles and add the edges to the egoGraph
            egoGraph.forNodes([&](const node localV) {
                const node globalV = egoMapping.toGlobal(localV);
                directedG.forEdgesOf(globalV, [&](node, node w, edgeweight weight) {
                    if (egoMapping.isMapped(w)) {
                        // we have found a triangle u-v-w
                        egoGraph.addEdge(localV, egoMapping.toLocal(w), weight);
                    }
                });
            });


            Partition egoPartition;

            {
                std::unique_ptr<CommunityDetectionAlgorithm> currentLocalClusteringAlgo(
                    localClusteringAlgo->clone());
                currentLocalClusteringAlgo->setGraph(egoGraph);
                currentLocalClusteringAlgo->run();
                egoPartition = currentLocalClusteringAlgo->getPartition();
            }

            if (egoPartition.upperBound() > egoGraph.numberOfNodes()) { // compact if necessary
                egoPartition.compact(egoPartition.upperBound() > 2 * egoGraph.numberOfNodes());
            }

            // Map parts that contain nodes of the original egoNet to persona indexes
            std::vector<index> partToPersona(egoPartition.upperBound(), none);
            const index egoPartsOffset = egoNetPartitionsOffset[egoNode];
            index numPersonas = 0;
            for (node u = 0; u < originalEgoNetSize; ++u) {
                index p = egoPartition[u];
                if (partToPersona[p] == none) {
                    partToPersona[p] = numPersonas;
                    ++numPersonas;
                }

                const node globalU = egoMapping.toGlobal(u);
                egoNetPartitions[egoPartsOffset + u] = std::make_pair(globalU, partToPersona[p]);
            }
            std::sort(egoNetPartitions.begin() + egoPartsOffset,
                      egoNetPartitions.begin() + egoPartsOffset + originalEgoNetSize);
            egoNetPartitionCounts[egoNode] = numPersonas;

            if (storeEgoNetGraphs) { // only for analysis
                storeEgoNetGraph(egoGraph, egoMapping, egoNode, egoPartition);
            }

            if (connectPersonas) {
                ParallelPartitionCoarsening coarsening(egoGraph, egoPartition, false);
                coarsening.run();
                const Graph &coarseGraph = coarsening.getCoarseGraph();
                const std::vector<node> &egoToCoarse = coarsening.getFineToCoarseNodeMapping();

                // Build a mapping from nodes in the coarse graph to partition ids in egoPartition
                std::vector<node> personaIndex(coarseGraph.upperNodeIdBound(), none);
                for (node u = 0; u < originalEgoNetSize; ++u) {
                    personaIndex[egoToCoarse[u]] = partToPersona[egoPartition[u]];
                    assert(personaIndex[egoToCoarse[u]] != none);
                }

                // Insert edges between the personas
                RandomMaximumSpanningForest span{coarseGraph};
                span.run();
                const Graph &spanningForest = span.getMSF();
                std::vector<WeightedEdge> edges;
                edges.reserve(std::min(spanningForest.numberOfEdges(), numPersonas - 1));
                spanningForest.forEdges([&](node u, node v) {
                    index p_u = personaIndex[u];
                    index p_v = personaIndex[v];
                    if (p_u != p_v && p_u != none && p_v != none)
                        edges.emplace_back(p_u, p_v, 1);
                });

                personaEdges[egoNode] = std::move(edges);
            }
        }
    }
}

void
EgoSplitting::storeEgoNetGraph(const Graph &egoGraph, const NodeMapping &egoMapping, node egoNode,
                               const Partition &egoPartition) {
    // Only store a given maximum of ego-nets (expected)
    count maxEgoNets = std::stoi(parameters.at("Max Stored EgoNets"));
    double storeChance = (double) maxEgoNets / G->numberOfNodes();
    if (Aux::Random::real() > storeChance)
        return;

    DEBUG("Store ego-net of node ", egoNode);
    // Get EgoNet with global node ids
    std::vector<WeightedEdge> edges;
    edges.reserve(egoGraph.numberOfEdges());
    // Add loops to retain isolated nodes
    egoGraph.forNodes([&](node u) {
        node globalId = egoMapping.toGlobal(u);
        edges.emplace_back(globalId, globalId, defaultEdgeWeight);
    });
    egoGraph.forEdges([&](node u, node v, edgeweight weight) {
        edges.emplace_back(egoMapping.toGlobal(u), egoMapping.toGlobal(v), weight);
    });
    egoNets[egoNode] = edges;

    // Store ego-net partition with extended nodes
    egoGraph.forNodes([&](node i) {
        egoNetExtendedPartitions[egoNode].emplace(egoMapping.toGlobal(i),
                                                  egoPartition.subsetOf(i));
    });
}

void EgoSplitting::splitIntoPersonas() {
    count sum = 0;
    G->forNodes([&](node u) {
        personaOffsets[u] = sum;
        sum += egoNetPartitionCounts[u];
    });
    personaOffsets[G->upperNodeIdBound()] = sum;
    personaGraph = Graph(sum, G->isWeighted());
}

void EgoSplitting::connectPersonas() {
    auto getPersona = [&](node u, index i) -> node {
        assert(i < egoNetPartitionCounts[u]);
        return personaOffsets[u] + i;
    };

    Aux::Timer timer;
    timer.start();

    for (node u = 0; u < G->upperNodeIdBound(); ++u) {
        if (!G->hasNode(u) || G->degree(u) < 2)
            continue;
        // Connect personas of each node
        for (const WeightedEdge &edge : personaEdges[u]) {
            node pu = getPersona(u, edge.u);
            node pv = getPersona(u, edge.v);
            assert(pu != pv);
            personaGraph.addEdge(pu, pv, edge.weight);
        }

        // Connect personas of different nodes
        const index egoBeginU = egoNetPartitionsOffset[u];
        const index nextEgoBegin = egoNetPartitionsOffset[u + 1];
        const index personaOffsetU = getPersona(u, 0);
        for (index i = egoBeginU; i < nextEgoBegin; ++i) {
            const std::pair<node, index> &neighborPersona = egoNetPartitions[i];
            node neighbor = neighborPersona.first;

            // Are we already at the end of the current ego node
            if (neighbor <= u || neighbor == none) break;

            index neighborPos = egoNetPartitionsOffset[neighbor];
            assert(neighborPos < egoNetPartitions.size());
            assert(egoNetPartitions[neighborPos].first == u);
            const node pu = personaOffsetU + neighborPersona.second;
            const node pv = getPersona(neighbor, egoNetPartitions[neighborPos].second);
            assert(pu != none && pv != none);
            personaGraph.addEdge(pu, pv);

            ++egoNetPartitionsOffset[neighbor];
            ++egoNetPartitionsOffset[u];
        }
    }

    timer.stop();
    INFO("adding persona graph edges took ", timer.elapsedMilliseconds(), " ms");

#ifndef NDEBUG
    {
        count internalPersonaEdges = 0;
        for (const auto &edges : personaEdges)
            internalPersonaEdges += edges.size();
        count removedEdges = 0;
        G->forNodes([&](node v) {
            if (G->degree(v) == 1) {
                node neighbor = none;
                G->forNeighborsOf(v, [&](node x) { neighbor = x; });
                assert(neighbor != none);
                if (G->degree(neighbor) > 1 || v < neighbor) {
                    ++removedEdges;
                }
            }
        });
        assert(personaGraph.numberOfEdges() + removedEdges == G->numberOfEdges() + internalPersonaEdges);

        personaGraph.forNodes([&](node u) {
            assert(personaGraph.degree(u) >= 1);
            if (personaGraph.degree(u) == 1) {
                node neighbor = none;
                personaGraph.forNeighborsOf(u, [&](node x) { neighbor = x; });
                assert(neighbor != u);
            }
        });

        // check that no isolated nodes were added
        ConnectedComponents compsAlgo(personaGraph);
        compsAlgo.run();
        auto componentSizeMap = compsAlgo.getPartition().subsetSizeMap();
        count isolatedNodes = 0;
        for (const auto &component  : componentSizeMap) {
            count numNodes = component.second;
            if (numNodes == 1)
                ++isolatedNodes;
        }
        assert(isolatedNodes == 0);
    }
#endif

    timer.start();

    // Remove small components in the persona graph, as the resulting communities would always be discarded
    ConnectedComponents compsAlgo(personaGraph);
    compsAlgo.run();
    const Partition& comps = compsAlgo.getPartition();

    timer.stop();
    INFO("connected components took ", timer.elapsedMilliseconds(), " ms");
    timer.start();

    std::vector<count> componentSizes(compsAlgo.numberOfComponents());
    personaGraph.parallelForNodes([&](node u) {
        if (componentSizes[comps[u]] < minCommunitySize) {
#pragma omp atomic update
            componentSizes[comps[u]] += 1;
        }

    });

    timer.stop();
    INFO("component sizes took ", timer.elapsedMilliseconds(), " ms");
    timer.start();

    count removedEdges = 0;
    personaGraph.forNodes([&](node u) {
        // the other direction gets deleted as well, since we remove connected components
        removedEdges += personaGraph.degree(u);
        personaGraph.removePartialOutEdges(unsafe, u);
        personaGraph.removeNode(u);
    });
    personaGraph.setEdgeCount(unsafe, personaGraph.numberOfEdges() - removedEdges / 2);

    timer.stop();
    INFO("removing small components took ", timer.elapsedMilliseconds(), " ms");

    egoNetPartitions.clear();
    egoNetPartitions.shrink_to_fit();
}

void EgoSplitting::createPersonaClustering() {
    globalClusteringAlgo->setGraph(personaGraph);
    globalClusteringAlgo->run();
    personaPartition = globalClusteringAlgo->getPartition();
}

std::vector<std::vector<node>> EgoSplitting::getCommunitiesFromPersonaClustering() {
    personaPartition.compact(personaPartition.upperBound() < 2 * personaGraph.numberOfNodes());
    std::vector<std::vector<node>> communities(personaPartition.upperBound());
    std::vector<std::atomic<index>> offsets(personaPartition.upperBound() + 1);
    personaGraph.balancedParallelForNodes([&](node persona) {
        offsets[personaPartition.subsetOf(persona)].fetch_add(1, std::memory_order_relaxed);
    });

    // no need to parallelize. it's already fast
    for (index i = 1; i < offsets.size(); ++i) {
        offsets[i] += offsets[i - 1];
    }

    // can actually use number of nodes here, not upperNodeIdBound !
    std::vector<node> sorted(personaGraph.numberOfNodes(), none);
    G->balancedParallelForNodes([&](node u) {
        for (index i = personaOffsets[u]; i < personaOffsets[u + 1]; ++i) {
            if (personaGraph.hasNode(i)) {
                index comm = personaPartition.subsetOf(i);
                index pos = offsets[comm].fetch_sub(1, std::memory_order_relaxed) - 1;
                sorted[pos] = u;
            }
        }
    });

#pragma omp parallel for schedule (guided)
    for (omp_index i = 0; i < offsets.size() - 1; ++i) {
        // each community sub-range may now have duplicate nodes!
        auto begin = sorted.begin() + offsets[i].load(std::memory_order_relaxed);
        auto end = sorted.begin() + offsets[i + 1].load(std::memory_order_relaxed);
        std::sort(begin, end);
        auto newEnd = std::unique(begin, end);
        auto& nodeList = communities[i];
        for ( ; begin != newEnd; ++begin) {
            nodeList.push_back(*begin);
        }
    }


#ifndef NDEBUG
    // check if each node in the node list of a community has a persona in that community
#pragma omp parallel for schedule (guided)
    for (omp_index c = 0; c < communities.size(); ++c) {
        for (node u : communities[c]) {
            bool any = false;
            for (node i = personaOffsets[u]; !any && i < personaOffsets[u + 1]; ++i) {
                any |= personaGraph.hasNode(i) && personaPartition.subsetOf(i) == c;
            }
            assert(any);
        }
    }

    // check if each node is in the node-lists of the communities of its personas
    G->balancedParallelForNodes([&](node u) {
      for (node i = personaOffsets[u]; i < personaOffsets[u + 1]; ++i) {
          if (personaGraph.hasNode(i)) {
              const index c = personaPartition.subsetOf(i);
              assert(std::find(communities[c].begin(), communities[c].end(), u) != communities[c].end());
          }
      }
    });
#endif
    return communities;
}

void EgoSplitting::createCover(const std::vector<std::vector<node>> &communities) {
    resultCover = Cover(G->upperNodeIdBound());
    resultCover.setUpperBound(communities.size());
    for (index com_id = 0; com_id < communities.size(); ++com_id) {
        for (node u : communities[com_id]) {
            resultCover.addToSubset(com_id, u);
        }
    }
}

void EgoSplitting::cleanUpCommunities(std::vector<std::vector<node>> &communities) {
    if (parameters.at("Cleanup") == "Yes") {
        discardSmallCommunities(communities);

        if (parameters.at("Cleanup Conductance") == "Yes") {
            // Cleanup does initialization potentially in O(m) - initialize once, copy O(n)
            // information for each thread.
            ConductanceCommunityCleanup cleanup(*G);

#pragma omp parallel for schedule(dynamic, 10) firstprivate(cleanup)
            for (omp_index communityIndex = 0; communityIndex < communities.size();
                 ++communityIndex) {
                std::vector<node> &community(communities[communityIndex]);
                TRACE("Starting cleanup of community of size ", community.size());

                cleanup.setCommunity(community.begin(), community.end());
                cleanup.run();

                community.clear();
                cleanup.forCommunityNodes([&](node u) { community.push_back(u); });
                assert(community.size() == cleanup.size());
            }
        } else {
            bool mergeDiscarded = parameters.at("Cleanup Merge") == "Yes";
            double minOverlap = Aux::stringToDouble(parameters.at("Cleanup Min Overlap"));
            SignificanceCommunityCleanUp cleanup(*G, communities, stochasticDistribution, 0.1, 0.1,
                                                 minOverlap, mergeDiscarded);
            cleanup.run();
        }
    }
    discardSmallCommunities(communities);

    INFO("Got ", communities.size(), " communities after cleanup");
}

void EgoSplitting::discardSmallCommunities(std::vector<std::vector<node>> &communities) {
    auto new_end = std::remove_if(communities.begin(), communities.end(),
                                  [&](const std::vector<node> &c) {
                                      return c.size() < minCommunitySize;
                                  });
    communities.erase(new_end, communities.end());
}

Cover EgoSplitting::getCover() {
    return resultCover;
}

std::string EgoSplitting::toString() const {
    return "EgoSplitting";
}

std::unordered_map<node, std::vector<WeightedEdge>> EgoSplitting::getEgoNets() {
    return egoNets;
}

std::vector<std::unordered_map<node, index>> EgoSplitting::getEgoNetPartitions() {
    return egoNetExtendedPartitions;
}

void
EgoSplitting::setParameters(std::map<std::string, std::string> const &new_parameters) {
    for (auto &x : new_parameters) {
        parameters.at(x.first) = x.second;
    }
}

} /* namespace NetworKit */
