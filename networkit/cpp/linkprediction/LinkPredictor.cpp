/*
 * LinkPredictor.cpp
 *
 *  Created on: 28.02.2015
 *      Author: Kolja Esders (kolja.esders@student.kit.edu)
 */

#include <algorithm>

#include <networkit/linkprediction/LinkPredictor.hpp>
#include <networkit/auxiliary/Log.hpp>
#include <networkit/auxiliary/Parallel.hpp>

#include <omp.h>

namespace NetworKit {

LinkPredictor::LinkPredictor() : G(nullptr), validCache(false) {
}

LinkPredictor::LinkPredictor(const Graph& G) : G(&G), validCache(false) {
}

void LinkPredictor::setGraph(const Graph& newGraph) {
  G = &newGraph;
  validCache = false;
}

std::vector<LinkPredictor::prediction> LinkPredictor::runOn(std::vector<std::pair<node, node>> nodePairs) {
  std::vector<prediction> predictions(nodePairs.size());
  Aux::Parallel::sort(nodePairs.begin(), nodePairs.end());
  #pragma omp parallel for schedule(dynamic) shared(predictions)
  for (omp_index i = 0; i < static_cast<omp_index>(nodePairs.size()); ++i) {
    predictions[i] = std::make_pair(nodePairs[i], run(nodePairs[i].first, nodePairs[i].second));
  }
  return predictions;
}

double LinkPredictor::run(node u, node v) {
  if (G == nullptr) {
    throw std::logic_error("Set a graph first.");
  } else if (!G->hasNode(u) || !G->hasNode(v)) {
    throw std::invalid_argument("Invalid node provided.");
  } else if (G->isDirected()) {
    throw std::invalid_argument("Only undirected graphs accepted.");
  } else if (u == v) {
    // Default behaviour for self-loops
    return 0;
  }
  return runImpl(u, v);
}

std::vector<LinkPredictor::prediction> LinkPredictor::runAll() {
  std::vector<node> nodes = G->nodes();
  std::vector<std::pair<node, node>> nodePairs;
  // Exclude all node-pairs that are already connected and ensure u != v for all node-pairs (u, v).
  for (index i = 0; i < nodes.size(); ++i) {
    for (index j = i + 1; j < nodes.size(); ++j) {
      if (!G->hasEdge(i, j))
        nodePairs.push_back(std::make_pair(i, j));
    }
  }
  return runOn(nodePairs);
}

} // namespace NetworKit