#include "CKBDynamic.h"


namespace NetworKit {

  namespace {
    std::pair<node, node> canonicalEdge(node u, node v) {
      if (u < v) return std::make_pair(u, v);
      return std::make_pair(v, u);
    }

    std::pair<node, node> edgeFromIndex(index i) {
      node u = 1 + std::floor(-0.5 + std::sqrt(0.25 + 2.0 * i));
      node v = i - (u * (u - 1) / 2);

      return canonicalEdge(u, v);
    }

    /**
     * Returns number of steps you need to wait until the next success (edge) occurs.
     */
    count get_next_edge_distance(const double log_cp) {
      return static_cast<count>(1 + floor(log(1.0 - Aux::Random::probability()) / log_cp));
    }

    class UniqueSampler {
    public:
      UniqueSampler(count max) : dist(0, max-1), max(max), i(0) {};
      count draw() {
	dist.param(std::uniform_int_distribution<count>::param_type(i, max-1));
	const count r = dist(Aux::Random::getURNG());
	count result = r;

	auto rit = replace.find(r);
	if (rit != replace.end()) {
	  result = rit->second;
	}

	auto iit = replace.find(i);
	if (iit == replace.end()) {
	  replace[r] = i;
	} else {
	  replace[r] = iit->second;
	}

	++i;

	return result;
      };
    private:
      std::uniform_int_distribution<count> dist;
      count max, i;
      tlx::btree_map<count, count> replace;
    };
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

  count CKBDynamic::Community::drawDesiredNumberOfEdges(double prob) const {
      std::binomial_distribution<count> distr(getMaximumNumberOfEdges(), prob);
      return distr(Aux::Random::getURNG());
  }

  void CKBDynamic::Community::perturbEdges(double prob) {
    assert(prob <= 1);
    assert(prob >= 0);

    const count desiredNumberOfEdges = drawDesiredNumberOfEdges(edgeProbability);
    std::binomial_distribution<count> distr(edges.size(), prob);
    const count edgesToPerturb = distr(Aux::Random::getURNG());
    count numEdgesToAdd = edgesToPerturb, numEdgesToRemove = edgesToPerturb;

    if (desiredNumberOfEdges > edges.size()) {
      numEdgesToAdd += desiredNumberOfEdges - edges.size();
    } else if (edges.size() < desiredNumberOfEdges) {
      numEdgesToRemove += edges.size() - desiredNumberOfEdges;
    }

    count numNonEdges = getMaximumNumberOfEdges() - edges.size();
    if (numEdgesToAdd > numNonEdges) {
      numEdgesToRemove -= (numEdgesToAdd - numNonEdges);
      numEdgesToAdd = numNonEdges;
    }

    std::vector<std::pair<node, node>> edgesToRemove, edgesToAdd;

    {
      UniqueSampler edgeSampler(edges.size());
      for (count i = 0; i < numEdgesToRemove; ++i) {
	auto e = edges.at(edgeSampler.draw());
	edgesToRemove.push_back(e);
      }
    }

    if (storeNonEdges) {
      UniqueSampler edgeSampler(nonEdges.size());
      for (count i = 0; i < numEdgesToAdd; ++i) {
	auto e = nonEdges.at(edgeSampler.draw());
	edgesToAdd.push_back(e);
      }
    } else {
      UniqueSampler edgeSampler(getMaximumNumberOfEdges());
      while (edgesToAdd.size() < numEdgesToAdd) {
	count edgeIndex = edgeSampler.draw();
	auto e = edgeFromIndex(edgeIndex);

	if (!edges.contains(e)) {
	  edgesToAdd.push_back(e);
	}
      }
    }

    for (auto e : edgesToAdd) {
      addEdge(e.first, e.second);
    }

    for (auto e : edgesToRemove) {
      removeEdge(e.first, e.second);
    }
  }

  void CKBDynamic::Community::removeRandomEdges(count k) {
    assert(k >= edges.size());

    for (count i = 0; i < k; ++k) {
      auto e = edges.at(Aux::Random::index(edges.size()));
      removeEdge(e.first, e.second);
    }
  }

  void CKBDynamic::Community::addRandomEdges(count k) {
    assert(edges.size() + k <= getMaximumNumberOfEdges());

    if (storeNonEdges) {
      assert(k >= nonEdges.size());

      for (count i = 0; i < k; ++k) {
	auto e = nonEdges.at(Aux::Random::index(nonEdges.size()));
	addEdge(e.first, e.second);
      }
    } else {
      const count numEdgesWanted = edges.size() + k;

      while (edges.size() < numEdgesWanted) {
	node u = nodes.at(Aux::Random::index(nodes.size()));
	node v = nodes.at(Aux::Random::index(nodes.size()));
	if (u != v && !edges.contains(canonicalEdge(u, v))) {
	  addEdge(u, v);
	}
      }
    }
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

    edgeProbability = prob;

    if (nodes.size() == 1) {
      return;
    }

    count desiredNumberOfEdges = drawDesiredNumberOfEdges(prob);

    if (desiredNumberOfEdges < edges.size()) {
      addRandomEdges(edges.size() - desiredNumberOfEdges);
    } else if (desiredNumberOfEdges > edges.size()) {
      removeRandomEdges(desiredNumberOfEdges - edges.size());
    }
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
