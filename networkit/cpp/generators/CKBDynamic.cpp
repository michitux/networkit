#include "CKBDynamic.h"


namespace NetworKit {

  namespace {
    std::pair<node, node> canonicalEdge(node u, node v) {
      if (u < v) return std::make_pair(u, v);
      return std::make_pair(v, u);
    }

    /**
     * Returns number of steps you need to wait until the next success (edge) occurs.
     */
    count get_next_edge_distance(const double log_cp) {
      return static_cast<count>(1 + floor(log(1.0 - Aux::Random::probability()) / log_cp));
    }
  }

  CKBDynamic::Community::Community(const Community& o) : edges(o.edges), nonEdges(o.nonEdges), nodes(o.nodes), neighbors(o.neighbors), edgeProbability(o.edgeProbability), isAvailable(o.isAvailable), storeNonEdges(o.storeNonEdges), generator(o.generator) {
    for (auto e : edges) {
      generator.addEdge(e.first, e.second);
    }

    for (auto u : nodes) {
      generator.addNodeToCommunity(u, CommunityPtr(this));
    }
  }


  CKBDynamic::Community::Community(double edgeProbability, CKBDynamic& generator) : edgeProbability(edgeProbability), isAvailable(false), storeNonEdges(edgeProbability > 0.6), generator(generator) {}

  void CKBDynamic::Community::removeEdge(node u, node v) {
      edges.erase(canonicalEdge(u, v));
      if (storeNonEdges) {
	nonEdges.insert(canonicalEdge(u, v));
      }
      neighbors[v].erase(u);
      neighbors[u].erase(v);
      generator.removeEdge(u, v);
  }

  void CKBDynamic::Community::addEdge(node u, node v) {
      neighbors[u].insert(v);
      neighbors[v].insert(u);
      edges.insert(canonicalEdge(u, v));
      if (storeNonEdges) {
	nonEdges.erase(canonicalEdge(u, v));
      }
      generator.addEdge(u, v);
  }

  void CKBDynamic::Community::removeNode(node u) {
    if (!nodes.contains(u)) throw std::runtime_error("Node not in community!");

    for (node v : neighbors[u]) {
      // We cannot delete from neighbors[u] because we iterate over it!
      edges.erase(canonicalEdge(u, v));
      neighbors[v].erase(u);
      generator.removeEdge(u, v);
    }

    neighbors.erase(u);
    nodes.erase(u);
    generator.removeNodeFromCommunity(u, CommunityPtr(this));

    if (storeNonEdges) {
      for (node i = 0; i < nodes.size(); ++i) {
	nonEdges.erase(canonicalEdge(u, nodes.at(i)));
      }
    }
  }

  void CKBDynamic::Community::addNode(node u) {
    if (nodes.contains(u)) throw std::runtime_error("Node already in community!");

    if (storeNonEdges) {
      for (node i = 0; i < nodes.size(); ++i) {
	nonEdges.insert(canonicalEdge(u, nodes.at(i)));
      }
    }

    neighbors.insert2(u, Aux::SamplingSet<node>());
    const double log_cp = std::log(1.0 - edgeProbability);

    for (node next = get_next_edge_distance(log_cp) - 1; next < nodes.size(); next += get_next_edge_distance(log_cp)) {
      addEdge(u, nodes.at(next));
    }

    nodes.insert(u);
  }

  void CKBDynamic::Community::changeEdgeProbability(double prob) {
    if (prob > 1) prob = 1;

    if (prob < 0.4 && storeNonEdges) {
      nonEdges.clear();
      storeNonEdges = false;
    } else if (prob > 0.6 && !storeNonEdges) {

      for (index i = 0; i < nodes.size(); ++i) {
	const node u = nodes.at(i);
	for (index j = i; j < nodes.size(); ++j) {
	  const node v = nodes.at(j);

	  const auto e = canonicalEdge(u, v);
	  if (!edges.contains(e)) {
	    nonEdges.insert(e);
	  }
	}
      }

      storeNonEdges = true;
    }

    if (nodes.size() == 1) {
      edgeProbability = prob;
      return;
    }

    const double actualEdgeProbability = edges.size() / (nodes.size() * (nodes.size() - 1));

    if (prob > actualEdgeProbability) {
      // insert edges
      if (storeNonEdges) {
	assert(nonEdges.size() + edges.size() == (nodes.size() * (nodes.size() - 1)));
	const double addProbability = 1.0 - ((1.0 - prob) / (1.0 - actualEdgeProbability));
	const double log_cp = std::log(1.0 - addProbability);

	std::vector<std::pair<node, node>> edges_to_add;

	for (index i = get_next_edge_distance(log_cp) - 1; i < nonEdges.size(); i += get_next_edge_distance(log_cp)) {
	  edges_to_add.push_back(nonEdges.at(i));
	}

	for (auto e : edges_to_add) {
	  addEdge(e.first, e.second);
	}
      } else {
	// insert edges without nonEdges
	// FIXME this should be randomized according to the G(n, p) model!
	const count numEdgesWanted = prob * nodes.size() * (nodes.size() - 1) / 2;

	while (edges.size() < numEdgesWanted) {
	  node u = nodes.at(Aux::Random::index(nodes.size()));
	  node v = nodes.at(Aux::Random::index(nodes.size()));
	  if (u != v && !edges.contains(canonicalEdge(u, v))) {
	    addEdge(u, v);
	  }
	}
      }
    } else if (prob < actualEdgeProbability) {
      // remove edges
      const double removeProbability = 1.0 - (prob / actualEdgeProbability);
      const double log_cp = std::log(1.0 - removeProbability);

      std::vector<std::pair<node, node>> edges_to_remove;

      for (index i = get_next_edge_distance(log_cp) - 1; i < edges.size(); i += get_next_edge_distance(log_cp)) {
	edges_to_remove.push_back(edges.at(i));
      }

      for (auto e : edges_to_remove) {
	removeEdge(e.first, e.second);
      }
    }

    edgeProbability = prob;
  }

  void CKBDynamic::Community::combineWith(Community& other) {
    assert(other.nodes == nodes);

    for (auto e : other.edges) {
      if (edges.contains(e)) {
	generator.removeEdge(e.first, e.second);
      } else {
	edges.insert(e);
	neighbors[e.first].insert(e.second);
	neighbors[e.second].insert(e.first);
      }
    }

    other.edges.clear();
    other.neighbors.clear();

    for (node u : other.nodes) {
      generator.removeNodeFromCommunity(u, CommunityPtr(&other));
    }
  }
}
