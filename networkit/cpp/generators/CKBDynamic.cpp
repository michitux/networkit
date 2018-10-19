#include "CKBDynamic.h"
#include <tlx/unused.hpp>


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

  CKBDynamic::Community::Community(const Community& o) : id(o.generator.nextCommunityId()), edges(o.edges), nonEdges(o.nonEdges), nodes(o.nodes), neighbors(o.neighbors), edgeProbability(o.edgeProbability), isAvailable(o.isAvailable), storeNonEdges(o.storeNonEdges), generator(o.generator) {
    for (auto e : edges) {
      generator.addEdge(e.first, e.second);
    }

    for (auto u : nodes) {
      generator.addNodeToCommunity(u, CommunityPtr(this));
    }
  }


  CKBDynamic::Community::Community(double edgeProbability, CKBDynamic& generator) : id(generator.nextCommunityId()), edgeProbability(edgeProbability), isAvailable(false), storeNonEdges(edgeProbability > 0.6), generator(generator) {}

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

  void CKBDynamic::Community::removeRandomNode() {
    if (nodes.size() == 0) throw std::runtime_error("Error, no nodes in community!");

    node u = nodes.at(Aux::Random::index(nodes.size()));
    removeNode(u);
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
    generator.addNodeToCommunity(u, CommunityPtr(this));
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
    assert(other.nodes.size() == nodes.size());

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
    other.nonEdges.clear();

    for (node u : other.nodes) {
      generator.removeNodeFromCommunity(u, CommunityPtr(&other));
    }

    other.nodes.clear();
  }

  CKBDynamic::CommunityChangeEvent::CommunityChangeEvent(CKBDynamic& generator, count numSteps) : active(true), numSteps(numSteps), currentStep(0), generator(generator) {}

  void CKBDynamic::CommunityChangeEvent::adaptProbability(CommunityPtr com, double targetProb) {
    double prob = (com->getEdgeProbability() * (numSteps - currentStep - 1) + targetProb) / (numSteps - currentStep);
    com->changeEdgeProbability(prob);
  }

  CKBDynamic::CommunityBirthEvent::CommunityBirthEvent(CommunityPtr community, std::vector<node> nodes, count coreSize, count numSteps, CKBDynamic& generator) : CommunityChangeEvent(generator, numSteps), nodes(std::move(nodes)), coreSize(coreSize), community(community) {
    std::shuffle(nodes.begin(), nodes.end(), Aux::Random::getURNG());
  }

  void CKBDynamic::CommunityBirthEvent::nextStep() {
    count numNodesToAdd = nodes.size() / (numSteps - currentStep);

    if (currentStep == 0) numNodesToAdd = coreSize;

    for (index i = 0; i < numNodesToAdd; ++i) {
      node u = nodes.back();
      nodes.pop_back();

      // skip nodes that have been removed in the meantime
      if (generator.hasNode(u)) {
	community->addNode(u);
      }
    }

    ++currentStep;

    if (currentStep == numSteps) {
      active = false;
      generator.addAvailableCommunity(community);
    }
  }

  CKBDynamic::CommunityDeathEvent::CommunityDeathEvent(CommunityPtr community, count coreSize, count numSteps, CKBDynamic& generator) : CommunityChangeEvent(generator, numSteps), community(community), coreSize(coreSize) {}


  void CKBDynamic::CommunityDeathEvent::nextStep() {
    count numNodesToRemove;
    // ensure that in the last step exactly coreSize nodes are removed
    if (currentStep < numSteps) {
      numNodesToRemove = (community->getNumberOfNodes() - coreSize) / (numSteps - 1 - currentStep);
    } else {
      numNodesToRemove = community->getNumberOfNodes();
    }

    for (index i = 0; i < numNodesToRemove; ++i) {
      community->removeRandomNode();
    }

    ++currentStep;

    if (currentStep == numSteps) {
      active = false;
      generator.removeCommunity(community);
    }
  }

  CKBDynamic::CommunityMergeEvent::CommunityMergeEvent(CommunityPtr communityA, CommunityPtr communityB, double targetEdgeProbability, count numSteps, CKBDynamic& generator) : CommunityChangeEvent(generator, numSteps), targetEdgeProbability(targetEdgeProbability), communityA(communityA), communityB(communityB) {
    count overlapSize = 0;

    for (node u : communityA->getNodes()) {
      if (communityB->hasNode(u)) {
	++overlapSize;
      } else {
	nodesToAddToB.push_back(u);
      }
    }

    std::shuffle(nodesToAddToB.begin(), nodesToAddToB.end(), Aux::Random::getURNG());

    for (node u : communityB->getNodes()) {
      if (!communityA->hasNode(u)) {
	nodesToAddToA.push_back(u);
      }
    }

    std::shuffle(nodesToAddToA.begin(), nodesToAddToA.end(), Aux::Random::getURNG());

    targetEdgeProbabilityPerCommunity = 1 - std::sqrt(1 - targetEdgeProbability);
  }

  void CKBDynamic::CommunityMergeEvent::nextStep() {
    auto addNodes = [&](CommunityPtr com, std::vector<node>& nodes) {
		      count numNodesToAdd = nodes.size() / (numSteps - currentStep);
		      for (index i = 0; i < numNodesToAdd; ++i) {
			node u = nodes.back();
			nodes.pop_back();

			// skip nodes that have been removed in the meantime
			if (generator.hasNode(u)) {
			  com->addNode(u);
			}
		      }
		    };

    addNodes(communityA, nodesToAddToA);
    addNodes(communityB, nodesToAddToB);

    adaptProbability(communityA, targetEdgeProbabilityPerCommunity);
    adaptProbability(communityB, targetEdgeProbabilityPerCommunity);

    ++currentStep;
    if (currentStep == numSteps) {
      active = false;
      communityA->combineWith(*communityB);
      generator.removeCommunity(communityB);
      generator.addAvailableCommunity(communityA);
    }
  }

  CKBDynamic::CommunitySplitEvent::CommunitySplitEvent(CommunityPtr community, count targetSizeA, double targetEdgeProbabilityA, count targetSizeB, double targetEdgeProbabilityB, count numSteps, CKBDynamic& generator) : CommunityChangeEvent(generator, numSteps), targetEdgeProbability({targetEdgeProbabilityA, targetEdgeProbabilityB}), communities({community, CommunityPtr(new Community(*community))}) {
    std::vector<node> nodes;
    nodes.reserve(communities[0]->getNumberOfNodes());
    for (node u : communities[0]->getNodes()) {
      nodes.push_back(u);
    }

    std::shuffle(nodes.begin(), nodes.end(), Aux::Random::getURNG());

    while (targetSizeA + targetSizeB < nodes.size()) {
      // sample a random set of nodes to remove from both
      node u = nodes.back();
      nodes.pop_back();
      nodesToRemove[0].push_back(u);
      nodesToRemove[1].push_back(u);
    }

    // Calculate how many nodes need to be removed from the two communities.
    // Due to the previous loop, both numbers can sum up to at most nodes.size().
    // If they sum up to less, we just keep some nodes in the overlap.
    const std::array<count, 2> numNodesToRemove {nodes.size() - targetSizeA, nodes.size() - targetSizeB};

    for (count com = 0; com < 2; ++com) {
      for (count i = 0; i < numNodesToRemove[com]; ++i) {
	node u = nodes.back();
	nodes.pop_back();
	nodesToRemove[com].push_back(u);
      }
    }

    // TODO: add second (new) community to generator somehow?
  }

  void CKBDynamic::CommunitySplitEvent::nextStep() {
    for (count com = 0; com < 2; ++com) {
      const count numNodesToRemove = nodesToRemove[com].size() / (numSteps - currentStep);
      for (count i = 0; i < numNodesToRemove; ++i) {
	const node u = nodesToRemove[com].back();
	nodesToRemove[com].pop_back();
	communities[com]->removeNode(u);
      }

      adaptProbability(communities[com], targetEdgeProbability[com]);
    }

    ++currentStep;
    if (currentStep == numSteps) {
      active = false;
      generator.addAvailableCommunity(communities[0]);
      generator.addAvailableCommunity(communities[1]);
    }
  }

  void CKBDynamic::addEdge(node u, node v) {
    auto e = canonicalEdge(u, v);
    auto it = edgesAlive.find(e);

    if (it == edgesAlive.end()) {
      edgesAlive.insert2(e, 1);
      graphEvents.emplace_back(GraphEvent::EDGE_ADDITION, e.first, e.second);
    } else {
      edgesAlive[e] += 1;
    }
  }

  void CKBDynamic::removeEdge(node u, node v) {
    auto e = canonicalEdge(u, v);
    auto it = edgesAlive.find(e);
    if (it == edgesAlive.end()) {
      throw std::runtime_error("Error, removing edge that does not exist");
    }

    if (it->second > 1) {
      edgesAlive[e] -= 1;
    } else {
      edgesAlive.erase(it);
      graphEvents.emplace_back(GraphEvent::EDGE_REMOVAL, e.first, e.second);
    }
  }

  void CKBDynamic::addNodeToCommunity(node u, CommunityPtr com) {
    if (com != globalCommunity) {
      nodeCommunities[u].insert(com);
      communityEvents.emplace_back(CommunityEvent::NODE_JOINS_COMMUNITY, u, com->getId());
      ++currentCommunityMemberships;
    }
  }

  void CKBDynamic::removeNodeFromCommunity(node u, CommunityPtr com) {
    if (com != globalCommunity) {
      nodeCommunities[u].erase(com);
      communityEvents.emplace_back(CommunityEvent::NODE_LEAVES_COMMUNITY, u, com->getId());
      --currentCommunityMemberships;
    }
  }

  void CKBDynamic::addAvailableCommunity(CommunityPtr com) {
    availableCommunities.insert(com);
  }

  void CKBDynamic::removeCommunity(CommunityPtr com) {
    availableCommunities.erase(com);
  }

  index CKBDynamic::nextCommunityId() {
    index result = maxCommunityId;
    ++maxCommunityId;
    return result;
  }

  CKBDynamic::CKBDynamic(count n, count minCommunitySize, count maxCommunitySize, double communitySizeExponent, double minSplitRatio, count minCommunityMembership, count maxCommunityMembership, double communityMembershipExponent, double eventProbability, double intraCommunityEdgeProbability, double intraCommunityEdgeExponent, double epsilon, count numTimesteps) : communityNodeSampler(0, minCommunityMembership, maxCommunityMembership, communityMembershipExponent), communitySizeSampler(new PowerlawCommunitySizeDistribution(minCommunitySize, maxCommunitySize, communitySizeExponent, intraCommunityEdgeProbability, intraCommunityEdgeExponent, minSplitRatio)), n(n), eventProbability(eventProbability), epsilon(epsilon), numTimesteps(numTimesteps) {
  }

  std::vector<GraphEvent> CKBDynamic::getGraphEvents() const {
    this->assureFinished();
    return graphEvents;
  }

  std::vector<CommunityEvent> CKBDynamic::getCommunityEvents() const {
    this->assureFinished();
    return communityEvents;
  }

  void CKBDynamic::generateNode() {
      node u = communityNodeSampler.addNode();
      for (node v = u; v < communityNodeSampler.getNumberOfNodes(); ++v) {
        nodesAlive.insert(v);
        nodeCommunities.emplace_back();
        globalCommunity->addNode(v);
      }
  }

  void CKBDynamic::run() {
    if (hasRun) throw std::runtime_error("Error, run has already been called");

    // initialization
    globalCommunity = CommunityPtr(new Community(epsilon, *this));

    for (node u = 0; u < n; ++u) {
      generateNode();
    }

    while (currentCommunityMemberships < communityNodeSampler.getSumOfDesiredMemberships()) {
      count communitySize;
      double edgeProbability;
      std::tie(communitySize, edgeProbability) = communitySizeSampler->drawCommunity();
      std::vector<node> comNodes(communityNodeSampler.birthCommunityNodes(communitySize));
      CommunityPtr com(new Community(edgeProbability, *this));
      for (node u : comNodes) {
        com->addNode(u);
      }
      addAvailableCommunity(com);
    }

    for (count timestep = 0; timestep < numTimesteps; ++timestep) {
      graphEvents.emplace_back(GraphEvent::TIME_STEP);
      communityEvents.emplace_back(CommunityEvent::TIME_STEP);
    }

    graphEvents.emplace_back(GraphEvent::TIME_STEP);
    communityEvents.emplace_back(CommunityEvent::TIME_STEP);

    hasRun = true;
  }


  CKBDynamic::PowerlawCommunitySizeDistribution::PowerlawCommunitySizeDistribution(count minSize, count maxSize, double gamma, double alpha, double densityGamma, double minSplitRatio) : sequence(minSize, maxSize, gamma), alpha(alpha), densityGamma(densityGamma), minSplitRatio(minSplitRatio) {
    sequence.run();
  }

  double CKBDynamic::PowerlawCommunitySizeDistribution::getProbability(count size) const {
    double prob = alpha / std::pow(size, densityGamma);
    if (prob > 1 || size <= 2) prob = 1;
    return prob;
  }

  std::pair<count, double> CKBDynamic::PowerlawCommunitySizeDistribution::drawCommunity() {
		count size = sequence.getDegree();
		double probability = getProbability(size);
		return std::make_pair(size, probability);
	}

  std::pair<count, double> CKBDynamic::PowerlawCommunitySizeDistribution::mergeCommunities(count sizeA, double probabilityA, count sizeB, double probabilityB, count combinedNodes) {
    tlx::unused(sizeA);
    tlx::unused(probabilityA);
    tlx::unused(sizeB);
    tlx::unused(probabilityB);
		assert(getProbability(sizeA) == probabilityA);
		assert(getProbability(sizeB) == probabilityB);

		return std::make_pair(combinedNodes, getProbability(combinedNodes));
	}

  std::pair<std::pair<count, double>, std::pair<count, double>> CKBDynamic::PowerlawCommunitySizeDistribution::splitCommunity(count size, double) {
		count minSize = std::max<count>(minSplitRatio * size, sequence.getMinimumDegree());
		if (minSize * 2 > size) {
			throw std::runtime_error("Given community too small to split");
		}

		std::mt19937 gen(rand());
		std::uniform_int_distribution<> dis(minSize, size - minSize);
		int sizeA = dis(gen);
		int sizeB = size-sizeA;

		return std::make_pair(std::make_pair(sizeA, getProbability(sizeA)),std::make_pair(sizeB, getProbability(sizeB)));
	}

}
