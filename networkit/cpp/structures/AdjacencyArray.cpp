/*
 * AdjacencyArray.cpp
 *
 * Created: 2018-12-12
 * Author: Armin Wiebigke
 */

#include "AdjacencyArray.h"

namespace NetworKit {

AdjacencyArray::AdjacencyArray(const NetworKit::Graph &G) {
    edgesBegin.resize(G.numberOfNodes() + 1);
    edges.resize(G.numberOfEdges());

    // direct edge from low to high-degree nodes
    auto isOutEdge = [&](node u, node v) {
        return G.degree(u) < G.degree(v) || (G.degree(u) == G.degree(v) && u < v);
    };

    index pos = 0;
    for (index u = 0; u < G.upperNodeIdBound(); ++u) {
        edgesBegin[u] = pos;
        if (G.hasNode(u)) {
            G.forEdgesOf(u, [&](node, node v) {
                if (isOutEdge(u, v)) {
                    edges[pos++] = v;
                }
            });
        }
    }

    edgesBegin[G.upperNodeIdBound()] = pos;
}

} /* namespace NetworKit */