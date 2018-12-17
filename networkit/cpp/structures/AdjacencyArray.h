/*
 * AdjacencyArray.h
 *
 * Created: 2018-12-12
 * Author: Armin Wiebigke
 */

#ifndef ADJACENCYARRAY_H
#define ADJACENCYARRAY_H

#include "../graph/Graph.h"

namespace NetworKit {

// TODO: other name
class AdjacencyArray {
private:
    std::vector <index> edgesBegin;
    std::vector <node> edges;

public:
    AdjacencyArray(const Graph &G);

    template <typename L> void forEdgesOf(node u, L handle) const;

};


template<typename L>
void AdjacencyArray::forEdgesOf(node u, L handle) const {
    for (index i = edgesBegin[u]; i < edgesBegin[u + 1]; ++i) {
        handle(u, edges[i]);
    }
}

} /* namespace NetworKit */

#endif //ADJACENCYARRAY_H
