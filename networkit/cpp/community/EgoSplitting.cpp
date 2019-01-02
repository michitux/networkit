/*
 * EgoSplitting.h
 *
 * Created: 2018-12-11
 * Author: Armin Wiebigke
 */

#include <omp.h>
#include <iostream>
#include <cmath>

#include "EgoSplitting.h"
#include "../structures/UnionFind.h"
#include "../structures/Partition.h"
#include "../structures/AdjacencyArray.h"
#include "../components/ConnectedComponents.h"

namespace NetworKit {

EgoSplitting::EgoSplitting(const Graph &G,
                           const std::function<Partition(Graph &)> &localClusterAlgo,
                           const std::function<Partition(Graph &)> &globalClusterAlgo)
        : G(G), localClusterAlgo(localClusterAlgo), globalClusterAlgo(globalClusterAlgo) {
    egoNets.resize(G.upperNodeIdBound());
    personaOffsets.resize(G.upperNodeIdBound() + 1, 0);
}

void EgoSplitting::run() {
    createEgoNets();
    splitIntoPersonas();
    connectPersonas();
    createPersonaClustering();
    createCover();
}

std::string EgoSplitting::toString() const {
    return "EgoSplitting";
}

void EgoSplitting::createEgoNets() {
    AdjacencyArray directedEdges(G); // store each undirected edge as one directed edge
    // Assign IDs to the neighbours
    std::vector<std::vector<count> > nodeToId(omp_get_max_threads(),
                                              std::vector<count>(G.upperNodeIdBound(), none));

    // Store number of partitions of the ego-net
    for (std::size_t i = 0; i < G.upperNodeIdBound(); ++i) {
        egoNets[i].emplace(none, 0);
    }

    G.balancedParallelForNodes([&](node u) {
        auto tid = omp_get_thread_num();

        // Assign IDs from 0 to degree-1 to neighbors
        std::vector<node> idToNode(G.degree(u));
        {
            index i = 0;
            G.forEdgesOf(u, [&](node, node v) {
                idToNode[i] = v;
                nodeToId[tid][v] = i++;
            });
        }

        // Find all triangles and add the edges to the egoGraph
        Graph egoGraph(G.degree(u));
        G.forEdgesOf(u, [&](node, node v) {
            directedEdges.forEdgesOf(v, [&](node, node w) {
                if (nodeToId[tid][w] != none) {
                    // we have found a triangle u-v-w
                    egoGraph.addEdge(nodeToId[tid][v], nodeToId[tid][w]);
                }
            });
        });

        // Cluster ego-net with the local cluster algorithm
        auto egoPartition = localClusterAlgo(egoGraph);
        egoPartition.compact();

        // Insert nodes into ego-net data structure
        for (index i = 0; i < G.degree(u); ++i) {
            egoNets[u].emplace(idToNode[i], egoPartition.subsetOf(i));
        }
        egoNets[u][none] = egoPartition.numberOfSubsets();

        // Reset IDs
        for (node v : idToNode) {
            nodeToId[tid][v] = none;
        }
    });
}

void EgoSplitting::splitIntoPersonas() {
    count sum = 0;
    for (index i = 0; i < G.upperNodeIdBound(); ++i) {
        personaOffsets[i] = sum;
        //sum += std::max((count) 1, egoNets[i][none]); // include isolated nodes in persona graph?
        sum += egoNets[i][none];
    }
    personaOffsets[G.upperNodeIdBound()] = sum;
    personaGraph = Graph(sum);
}

void EgoSplitting::connectPersonas() {
    auto getPersona = [&](node u, index i) {
        return personaOffsets[u] + i;
    };

    G.forEdges([&](node u, node v) {
        auto idx_u = egoNets[u].find(v);
        auto idx_v = egoNets[v].find(u);
        assert(idx_u != egoNets[u].end() && idx_v != egoNets[v].end());
        personaGraph.addEdge(getPersona(u, idx_u->second), getPersona(v, idx_v->second));
    });

    egoNets.clear();

#ifndef NDEBUG
    assert(personaGraph.numberOfEdges() == G.numberOfEdges());
    // check that no isolated nodes were added
    auto numIsolatedNodes = [](const Graph &G) {
        ConnectedComponents compsAlgo(G);
        compsAlgo.run();
        auto comps = compsAlgo.getComponentSizes();
        count isolated = 0;
        for (const auto &x  : comps) {
            if (x.second == 1)
                ++isolated;
            assert(x.second != 0);
            assert(x.second != 1);
        }
        return isolated;
    };
    auto gEdges = G.edges();
    auto pEdges = personaGraph.edges();
    //auto iso1 = numIsolatedNodes(G);
    auto iso2 = numIsolatedNodes(personaGraph);
    assert(iso2 == 0);
#endif
}

void EgoSplitting::createPersonaClustering() {
    personaPartition = globalClusterAlgo(personaGraph);
    assert(personaPartition.upperBound() <= personaGraph.upperNodeIdBound());
}

void EgoSplitting::createCover() {
    cover = Cover(G.upperNodeIdBound());
    cover.setUpperBound(personaPartition.upperBound());
    G.forNodes([&](node u) {
        for (index i = personaOffsets[u]; i < personaOffsets[u + 1]; ++i) {
            cover.addToSubset(personaPartition.subsetOf(i), u);
        }
    });

#ifndef NDEBUG
    std::cout << "Detected Community Sizes:\n";
    for (auto s : cover.subsetSizeMap()) {
        std::cout << s.second << ", ";
    }
    std::cout << std::endl;
#endif

    // Discard communities of size 4 or less
    for (auto s : cover.subsetSizeMap()) {
        index set = s.first;
        count size = s.second;
        if (size <= 4) {
            cover.removeSubset(set);
        }
    }

#ifndef NDEBUG
    std::cout << "Result:\n";
    for (auto s : cover.subsetSizeMap()) {
        std::cout << s.second << ", ";
    }
    std::cout << std::endl;
#endif
}

Cover EgoSplitting::getCover() {
    return cover;
}

} /* namespace NetworKit */