/*
 * EgoSplitting.h
 *
 * Created: 2018-12-11
 * Author: Armin Wiebigke
 */

#include <omp.h>

#include "EgoSplitting.h"
#include "../structures/UnionFind.h"
#include "../structures/Partition.h"
#include "../structures/AdjacencyArray.h"
#include "PLP.h"

namespace NetworKit {

EgoSplitting::EgoSplitting(const NetworKit::Graph &G) : CommunityDetectionAlgorithm(G) {
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

Partition EgoSplitting::getPartition() {

}

std::string EgoSplitting::toString() const {

}

void EgoSplitting::createEgoNets() {
    AdjacencyArray directedEdges(G); // store each undirected edge as one directed edge
    // Assign IDs to the neighbours
    std::vector<std::vector<count> > nodeToId(omp_get_max_threads(),
                                              std::vector<count>(G.upperNodeIdBound(), none));

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

        // Stores ego-net partitions
        UnionFind unionFind(G.degree(u));

        // Find all triangles
        G.forEdgesOf(u, [&](node, node v) {
            directedEdges.forEdgesOf(v, [&](node, node w) {
                if (nodeToId[tid][w] != none) {
                    // we have found a triangle u-v-w
                    unionFind.merge(nodeToId[tid][v], nodeToId[tid][w]);
                }
            });
        });

        // Reset IDs
        for (const auto &v : idToNode) {
            nodeToId[tid][v] = none;
        }

        // Insert nodes into ego-net hash table
        auto partitions = unionFind.toPartition();
        for (index i = 0; i < G.degree(u); ++i) {
            egoNets[u].emplace(idToNode[i], partitions.subsetOf(i));
        }
    });
}

void EgoSplitting::splitIntoPersonas() {
    count sum = 0;
    for (index i = 0; i < G.upperNodeIdBound(); ++i) {
        personaOffsets[i] = sum;
        sum += egoNets[i].size();
    }
    personaOffsets[G.upperNodeIdBound()] = sum;
    personaGraph = Graph(sum);
}

void EgoSplitting::connectPersonas() {
    auto
    getPersona = [&](node u, index i) {
        return personaOffsets[u] + i;
    };

    G.forEdges([&](node u, node v) {
        auto idx_u = egoNets[u].find(v);
        auto idx_v = egoNets[v].find(u);
        if (idx_u != egoNets[u].end() && idx_v != egoNets[v].end()) {
            personaGraph.addEdge(getPersona(u, idx_u->second), getPersona(v, idx_v->second));
        }
    });
    assert(personaGraph.numberOfEdges() == G.numberOfEdges());
}

void EgoSplitting::createPersonaClustering() {
    PLP plp(personaGraph);
    plp.run();
    personaPartition = plp.getPartition();
    personaPartition.compact();
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
}

Cover EgoSplitting::getCover() {
    return cover;
}

} /* namespace NetworKit */