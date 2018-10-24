#include "CKBDynamic.h"
#include <tlx/unused.hpp>
#include "../auxiliary/UniqueSampler.h"
#include "../auxiliary/SignalHandling.h"


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

  CKBDynamic::Community::Community(const Community& o) : id(o.generator.nextCommunityId()), edges(o.edges), nonEdges(o.nonEdges), nodes(o.nodes), neighbors(o.neighbors), edgeProbability(o.edgeProbability), available(o.available), storeNonEdges(o.storeNonEdges), generator(o.generator) {
    for (auto e : edges) {
      generator.addEdge(e.first, e.second);
    }

    for (auto u : nodes) {
      generator.addNodeToCommunity(u, CommunityPtr(this));
    }

    generator.addCommunity(CommunityPtr(this));
  }


  CKBDynamic::Community::Community(double edgeProbability, CKBDynamic& generator) : id(generator.nextCommunityId()), edgeProbability(edgeProbability), available(false), storeNonEdges(edgeProbability > 0.6), generator(generator) {
    generator.addCommunity(CommunityPtr(this));
  }

  void CKBDynamic::Community::verifyInvariants() const {
    if (storeNonEdges) {
      for (auto e : edges) {
        assert(canonicalEdge(e.first, e.second) == e);
        assert(!nonEdges.contains(e));
        assert(nodes.contains(e.first));
        assert(nodes.contains(e.second));
      }

      for (auto e : nonEdges) {
        assert(canonicalEdge(e.first, e.second) == e);
        assert(nodes.contains(e.first));
        assert(nodes.contains(e.second));
      }

      assert(nonEdges.size() + edges.size() == getMaximumNumberOfEdges());
    }
  }

  void CKBDynamic::Community::removeEdge(node u, node v) {
      edges.erase(canonicalEdge(u, v));
      if (storeNonEdges) {
	nonEdges.insert(canonicalEdge(u, v));
        verifyInvariants();
      }
      neighbors[v].erase(u);
      neighbors[u].erase(v);
      generator.removeEdge(u, v);
  }

  void CKBDynamic::Community::addEdge(node u, node v) {
    assert(u != v);
    assert(nodes.contains(u));
    assert(nodes.contains(v));
      neighbors[u].insert(v);
      neighbors[v].insert(u);
      edges.insert(canonicalEdge(u, v));
      if (storeNonEdges) {
	nonEdges.erase(canonicalEdge(u, v));
        verifyInvariants();
      }
      generator.addEdge(u, v);
  }

  std::pair<node, node> CKBDynamic::Community::edgeFromIndex(index i) const {
    node u = 1 + std::floor(-0.5 + std::sqrt(0.25 + 2.0 * i));
    node v = i - (u * (u - 1) / 2);

    return canonicalEdge(nodes.at(u), nodes.at(v));
  }


  void CKBDynamic::Community::removeNode(node u) {
    assert(nodes.contains(u));
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
      verifyInvariants();
    }
  }

  void CKBDynamic::Community::removeRandomNode() {
    assert(nodes.size() > 0);
    if (nodes.size() == 0) throw std::runtime_error("Error, no nodes in community!");

    node u = nodes.at(Aux::Random::index(nodes.size()));
    removeNode(u);
  }

  void CKBDynamic::Community::addNode(node u) {
    assert(!nodes.contains(u));
    if (nodes.contains(u)) throw std::runtime_error("Node already in community!");
    nodes.insert(u);

    if (storeNonEdges) {
      for (node i = 0; i < nodes.size(); ++i) {
        const node v = nodes.at(i);
        if (u != v) {
          nonEdges.insert(canonicalEdge(u, v));
        }
      }

      verifyInvariants();
    }

    neighbors.insert2(u, Aux::SamplingSet<node>());
    const double log_cp = std::log(1.0 - edgeProbability);

    for (node next = get_next_edge_distance(log_cp) - 1; next < nodes.size(); next += get_next_edge_distance(log_cp)) {
      const node v = nodes.at(next);
      if (u != v) {
        addEdge(u, nodes.at(next));
      }
    }

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
      Aux::UniqueSampler edgeSampler(edges.size());
      for (count i = 0; i < numEdgesToRemove; ++i) {
	auto e = edges.at(edgeSampler.draw());
	edgesToRemove.push_back(e);
      }
    }

    if (storeNonEdges) {
      Aux::UniqueSampler edgeSampler(nonEdges.size());
      for (count i = 0; i < numEdgesToAdd; ++i) {
	auto e = nonEdges.at(edgeSampler.draw());
	edgesToAdd.push_back(e);
      }
    } else {
      Aux::UniqueSampler edgeSampler(getMaximumNumberOfEdges());
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

    verifyInvariants();
  }

  void CKBDynamic::Community::removeRandomEdges(count k) {
    assert(k <= edges.size());

    for (count i = 0; i < k; ++i) {
      auto e = edges.at(Aux::Random::index(edges.size()));
      removeEdge(e.first, e.second);
    }
  }

  void CKBDynamic::Community::addRandomEdges(count k) {
    assert(edges.size() + k <= getMaximumNumberOfEdges());

    if (storeNonEdges) {
      assert(k <= nonEdges.size());
      verifyInvariants();

      for (count i = 0; i < k; ++i) {
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
	for (index j = i + 1; j < nodes.size(); ++j) {
	  const node v = nodes.at(j);

	  const auto e = canonicalEdge(u, v);
	  if (!edges.contains(e)) {
	    nonEdges.insert(e);
	  }
	}
      }

      verifyInvariants();

      storeNonEdges = true;
    }

    edgeProbability = prob;

    if (nodes.size() == 1) {
      return;
    }

    count desiredNumberOfEdges = drawDesiredNumberOfEdges(prob);

    if (desiredNumberOfEdges < edges.size()) {
      removeRandomEdges(edges.size() - desiredNumberOfEdges);
    } else if (desiredNumberOfEdges > edges.size()) {
      addRandomEdges(desiredNumberOfEdges - edges.size());
    }

    assert(edges.size() == desiredNumberOfEdges);
  }

  void CKBDynamic::Community::combineWith(Community& other) {
    assert(&other != this);
    assert(other.nodes.size() == nodes.size());

    for (auto e : other.edges) {
      // First add edge here to ensure it exists at least once globally
      // so we don't generate remove/add events globally.
      if (!edges.contains(e)) {
        addEdge(e.first, e.second);
      }

      other.removeEdge(e.first, e.second);
    }

    other.edges.clear();
    other.neighbors.clear();
    other.nonEdges.clear();

    for (node u : other.nodes) {
      generator.removeNodeFromCommunity(u, CommunityPtr(&other));
    }

    other.nodes.clear();
  }

  void CKBDynamic::Community::setAvailable( bool avail ) {
    if (avail) {
      assert(nodes.size() > 0);
    }
    if (this->available != avail) {
      this->available = avail;
      generator.addCommunity(CommunityPtr(this));
    }
  }

  CKBDynamic::CommunityChangeEvent::CommunityChangeEvent(CKBDynamic& generator, count numSteps) : active(true), numSteps(numSteps), currentStep(0), generator(generator) {}

  void CKBDynamic::CommunityChangeEvent::adaptProbability(CommunityPtr com, double targetProb) {
    double prob = (com->getEdgeProbability() * (numSteps - currentStep - 1) + targetProb) / (numSteps - currentStep);
    com->changeEdgeProbability(prob);
  }

  CKBDynamic::CommunityBirthEvent::CommunityBirthEvent(count coreSize, count targetSize, double edgeProbability, count numSteps, CKBDynamic& generator) : CommunityChangeEvent(generator, numSteps), coreSize(coreSize), targetSize(targetSize), community(new Community(edgeProbability, generator)) {}

  CKBDynamic::CommunityBirthEvent::CommunityBirthEvent(count numSteps, CKBDynamic& generator) : CommunityChangeEvent(generator, numSteps)  {
    double edgeProbability;
    std::tie(targetSize, edgeProbability) = generator.communitySizeSampler->drawCommunity();
    coreSize = std::max<count>(0.1 * targetSize, generator.communitySizeSampler->getMinSize());
    community = CommunityPtr(new Community(edgeProbability, generator));
  }

  void CKBDynamic::CommunityBirthEvent::nextStep() {
    count numNodesToAdd = (targetSize - community->getNumberOfNodes()) / (numSteps - currentStep);

    // Ensure that after every step the community has at least coreSize nodes,
    // even if nodes have been deleted in the meantime.
    if (community->getNumberOfNodes() < coreSize) {
      numNodesToAdd = std::max(coreSize - community->getNumberOfNodes(), numNodesToAdd);
    }

    std::vector<node> nodesToAdd = generator.communityNodeSampler.birthCommunityNodes(numNodesToAdd, community->getNodes());

    // If we did not get enough nodes to make the community larger than the minimum size
    // let the community die.
    if (nodesToAdd.size() + community->getNumberOfNodes() < coreSize) {
      INFO("Sampler returned too few nodes, letting the community die.");
      active = false;
      while (community->getNumberOfNodes() > 0) {
        community->removeRandomNode();
      }
      generator.removeCommunity(community);
    } else {
      for (node u : nodesToAdd) {
        assert(generator.hasNode(u));
        community->addNode(u);
      }

      ++currentStep;

      if (currentStep == numSteps) {
        if (community->getNumberOfNodes() == targetSize) {
          active = false;
          community->setAvailable(true);
        } else {
          // Delay event completion because not enough nodes could be requested.
          INFO("Delaying community birth because not enough nodes could be sampled.");
          --currentStep;
        }
      }
    }
  }

  CKBDynamic::CommunityDeathEvent::CommunityDeathEvent(CommunityPtr community, count coreSize, count numSteps, CKBDynamic& generator) : CommunityChangeEvent(generator, numSteps), community(community), coreSize(coreSize) {
    community->setAvailable(false);
  }


  void CKBDynamic::CommunityDeathEvent::nextStep() {
    count numNodesToRemove;
    // Ensure that in the last step exactly coreSize nodes are removed.
    // If there are less than coreSize nodes remaining due to additional node
    // deletions let the community die earlier.
    if (currentStep < numSteps - 1 && community->getNumberOfNodes() >= coreSize) {
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

  CKBDynamic::CommunityMergeEvent::CommunityMergeEvent(CommunityPtr communityA, CommunityPtr communityB, count numSteps, CKBDynamic& generator) : CommunityChangeEvent(generator, numSteps), targetEdgeProbability(.0), communityA(communityA), communityB(communityB) {
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

    std::tie(std::ignore, targetEdgeProbability) = generator.communitySizeSampler->mergeCommunities(communityA->getNumberOfNodes(), communityA->getEdgeProbability(), communityB->getNumberOfEdges(), communityB->getEdgeProbability(), communityA->getNumberOfNodes() + communityB->getNumberOfNodes() - overlapSize);
    targetEdgeProbabilityPerCommunity = 1 - std::sqrt(1 - targetEdgeProbability);

    communityA->setAvailable(false);
    communityB->setAvailable(false);
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
      const count oldNodes = communityA->getNumberOfNodes();
      tlx::unused(oldNodes);
      communityA->combineWith(*communityB);
      generator.removeCommunity(communityB);
      assert(communityA->getNumberOfNodes() == oldNodes);
      // This shouldn't change much but otherwise the community will loose half of
      // the edges in the next perturbation.
      communityA->changeEdgeProbability(targetEdgeProbability);
      communityA->setAvailable(true);
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

    communities[0]->setAvailable(false);
    communities[1]->setAvailable(false);
  }

  void CKBDynamic::CommunitySplitEvent::nextStep() {
    for (count com = 0; com < 2; ++com) {
      const count numNodesToRemove = nodesToRemove[com].size() / (numSteps - currentStep);
      for (count i = 0; i < numNodesToRemove; ++i) {
	const node u = nodesToRemove[com].back();
	nodesToRemove[com].pop_back();
        if (generator.hasNode(u)) {
          communities[com]->removeNode(u);
        }
      }

      adaptProbability(communities[com], targetEdgeProbability[com]);
    }

    ++currentStep;
    if (currentStep == numSteps) {
      active = false;
      communities[0]->setAvailable(true);
      communities[1]->setAvailable(true);
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
      communityNodeSampler.assignCommunity(u);

      if (com->isAvailable()) {
        if (com->getNumberOfNodes() == 2*communitySizeSampler->getMinSize()) {
          splittableCommunities.insert(com);
        }

        if (com->getNumberOfNodes() > (communitySizeSampler->getMaxSize() - communitySizeSampler->getMinSize())) {
          mergeableCommunities.erase(com);
        }
      }
    }
  }

  void CKBDynamic::removeNodeFromCommunity(node u, CommunityPtr com) {
    if (com != globalCommunity) {
      nodeCommunities[u].erase(com);
      communityEvents.emplace_back(CommunityEvent::NODE_LEAVES_COMMUNITY, u, com->getId());
      communityNodeSampler.leaveCommunity(u);
      --currentCommunityMemberships;

      if (com->isAvailable()) {
        if (com->getNumberOfNodes() < 2*communitySizeSampler->getMinSize()) {
          splittableCommunities.erase(com);
        }

        if (com->getNumberOfNodes() <= communitySizeSampler->getMaxSize()) {
          mergeableCommunities.insert(com);
        }
      }
    }
  }

  void CKBDynamic::addCommunity(CommunityPtr com) {
    // FIXME: do not add global community? but impossible because addCommunity is called in constructor...
    if (com->isAvailable()) {
      // If community is too small, remove community again!!
      if (com->getNumberOfNodes() < communitySizeSampler->getMinSize()) {
        INFO("community has only ", com->getNumberOfNodes(), " nodes, destroying.");
        currentEvents.emplace_back(new CommunityDeathEvent(com, 0, 1, *this));
      } else {
        availableCommunities.insert(com);
        if (com->getNumberOfNodes() >= 2*communitySizeSampler->getMinSize()) {
          splittableCommunities.insert(com);
        }

        if (com->getNumberOfNodes() <= communitySizeSampler->getMaxSize() - communitySizeSampler->getMinSize()) {
          mergeableCommunities.insert(com);
        }
      }
    } else {
      availableCommunities.erase(com);
      splittableCommunities.erase(com);
      mergeableCommunities.erase(com);
    }
    communities.insert(com);
  }

  void CKBDynamic::removeCommunity(CommunityPtr com) {
    availableCommunities.erase(com);
    communities.erase(com);
    mergeableCommunities.erase(com);
    splittableCommunities.erase(com);
  }

  index CKBDynamic::nextCommunityId() {
    index result = maxCommunityId;
    ++maxCommunityId;
    return result;
  }

  CKBDynamic::CKBDynamic(count n, count minCommunitySize, count maxCommunitySize, double communitySizeExponent, double minSplitRatio, count minCommunityMembership, count maxCommunityMembership, double communityMembershipExponent, double communityEventProbability, double nodeEventProbability, double perturbationProbability, double intraCommunityEdgeProbability, double intraCommunityEdgeExponent, double epsilon, count numTimesteps) : communityNodeSampler(0, minCommunityMembership, maxCommunityMembership, communityMembershipExponent), communitySizeSampler(new PowerlawCommunitySizeDistribution(minCommunitySize, maxCommunitySize, communitySizeExponent, intraCommunityEdgeProbability, intraCommunityEdgeExponent, minSplitRatio)), n(n), communityEventProbability(communityEventProbability), nodeEventProbability(nodeEventProbability), perturbationProbability(perturbationProbability), epsilon(epsilon), numTimesteps(numTimesteps), currentCommunityMemberships(0) {
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

  void CKBDynamic::eraseNode() {
    node u = nodesAlive.at(Aux::Random::index(nodesAlive.size()));
    while (nodeCommunities[u].size() > 0) {
      CommunityPtr com = *nodeCommunities[u].begin();
      com->removeNode(u);
      // if a community becomes too small, erase it
      if (com->isAvailable() && com->getNumberOfNodes() < communitySizeSampler->getMinSize()) {
        INFO("Available community has only ", com->getNumberOfNodes(), " nodes, destroying.");
        currentEvents.emplace_back(new CommunityDeathEvent(com, 0, 1, *this));
      }
    }

    nodesAlive.erase(u);
    communityNodeSampler.removeNode(u);
    globalCommunity->removeNode(u);
  }

  count CKBDynamic::sampleNumSteps() const {
    return Aux::Random::integer(5, 15);
  }

  void CKBDynamic::run() {
    if (hasRun) throw std::runtime_error("Error, run has already been called");

    Aux::SignalHandler handler;

    // initialization
    globalCommunity = CommunityPtr(new Community(epsilon, *this));

    for (node u = 0; u < n; ++u) {
      generateNode();
    }

    const count initialNumberOfNodes = nodesAlive.size();

    while (currentCommunityMemberships < communityNodeSampler.getSumOfDesiredMemberships()) {
      handler.assureRunning();
      count communitySize;
      double edgeProbability;
      std::tie(communitySize, edgeProbability) = communitySizeSampler->drawCommunity();
      std::vector<node> comNodes(communityNodeSampler.birthCommunityNodes(communitySize));
      CommunityPtr com(new Community(edgeProbability, *this));
      for (node u : comNodes) {
        com->addNode(u);
      }
      com->setAvailable(true);
    }

    std::binomial_distribution<count> numEventDistribution;
    double deathProbability = 0.25, birthProbability = 0.25, splitProbability = 0.25, mergeProbability = 0.25;
    tlx::unused(mergeProbability);

    for (count timestep = 0; timestep < numTimesteps; ++timestep) {
      handler.assureRunning();
      graphEvents.emplace_back(GraphEvent::TIME_STEP);
      communityEvents.emplace_back(CommunityEvent::TIME_STEP);

      numEventDistribution.param(std::binomial_distribution<count>::param_type(communities.size(), communityEventProbability));
      const count numCommunityEvents = numEventDistribution(Aux::Random::getURNG());

      numEventDistribution.param(std::binomial_distribution<count>::param_type(nodesAlive.size(), nodeEventProbability));
      const count numNodeEvents = numEventDistribution(Aux::Random::getURNG());

      INFO("Timestep ", timestep, " generating ", numCommunityEvents, " community events and ", numNodeEvents, " node events");

      for (count i = 0; i < numCommunityEvents; ++i) {
        handler.assureRunning();
        count numSteps = sampleNumSteps();
        double r = Aux::Random::real();
        if (r < birthProbability) {
          // generate new community
          currentEvents.emplace_back(new CommunityBirthEvent(numSteps, *this));
        } else if (r < birthProbability + deathProbability) {
          // let a community die
          if (availableCommunities.size() > 0) {
            CommunityPtr com = availableCommunities.at(Aux::Random::index(availableCommunities.size()));
            count coreSize = std::max<count>(0.1 * com->getNumberOfNodes(), communitySizeSampler->getMinSize());
            currentEvents.emplace_back(new CommunityDeathEvent(com, coreSize, numSteps, *this));
            assert(!com->isAvailable());
          } else {
            WARN("No community available for death event.");
          }
        } else if (r < birthProbability + deathProbability + splitProbability) {
          // Split a community
          if (splittableCommunities.size() > 0) {
            CommunityPtr com = splittableCommunities.at(Aux::Random::index(splittableCommunities.size()));
            auto comSizeProb = communitySizeSampler->splitCommunity(com->getNumberOfNodes(), com->getEdgeProbability());
            currentEvents.emplace_back(new CommunitySplitEvent(com, comSizeProb.first.first, comSizeProb.first.second, comSizeProb.second.first, comSizeProb.second.second, numSteps, *this));
            assert(!com->isAvailable());
          } else {
            WARN("No community available for splitting.");
          }
        } else {
          // merge two communities
          if (mergeableCommunities.size() > 1) {
            CommunityPtr comA, comB;
            bool found = false;
            for (count j = 0; j < 20; ++j) {
              comA = mergeableCommunities.at(Aux::Random::index(mergeableCommunities.size()));
              comB = mergeableCommunities.at(Aux::Random::index(mergeableCommunities.size()));
              if (comA != comB && comA->getNumberOfNodes() + comB->getNumberOfNodes() < communitySizeSampler->getMaxSize()) {
                found = true;
                break;
              }
            }

            if (found) {
              currentEvents.emplace_back(new CommunityMergeEvent(comA, comB, numSteps, *this));
              assert(!comA->isAvailable());
              assert(!comB->isAvailable());
            } else {
              WARN("In 20 trials, no two communities found to merge.");
            }
          } else {
            WARN("No two communities available for merge.");
          }
        }
      } // generated all new community events

      // generate node events
      const double wantedNodeFraction = initialNumberOfNodes * 1.0 / nodesAlive.size();
      const double nodeBirthProbability = wantedNodeFraction / (1 + wantedNodeFraction);
      for (count j = 0; j < numNodeEvents; ++j) {
        if (Aux::Random::real() < nodeBirthProbability) {
          generateNode();
        } else {
          eraseNode();
        }
      }

      // Trigger all current events
      for (size_t e = 0; e < currentEvents.size();) {
        handler.assureRunning();
        currentEvents[e]->nextStep();

        if (!currentEvents[e]->isActive()) {
          std::swap(currentEvents[e], currentEvents.back());
          currentEvents.pop_back();
        } else {
          ++e;
        }
      }

      // edge perturbations
      if (perturbationProbability > 0) {
        globalCommunity->perturbEdges(perturbationProbability);

        const double sqrtPerturbationProbability = std::sqrt(perturbationProbability);

        const double log_perturb = std::log(1.0 - sqrtPerturbationProbability);

        for (count ci = get_next_edge_distance(log_perturb) - 1; ci < communities.size(); ci += get_next_edge_distance(log_perturb)) {
          handler.assureRunning();
          communities.at(ci)->perturbEdges(sqrtPerturbationProbability);
        }
      }

      // adjust event probabilities
      {
        const double x = communityNodeSampler.getSumOfDesiredMemberships() * 1.0 / currentCommunityMemberships;
        birthProbability = 0.5 * x / (1 + x);
        deathProbability = 0.5 - birthProbability;
        INFO("At timestep ", timestep, " adjusting birth probability to ", birthProbability, " and death probability to ", deathProbability);
        INFO("Current memberships: ", currentCommunityMemberships, " desired: ", communityNodeSampler.getSumOfDesiredMemberships(), " number of communities: ", communities.size(), " available: ", availableCommunities.size(), " active events ", currentEvents.size());
        INFO("Current nodes ", nodesAlive.size(), " current edges: ", edgesAlive.size(), " total graph events ", graphEvents.size(), " total community events ", communityEvents.size());
      }

    }

    graphEvents.emplace_back(GraphEvent::TIME_STEP);
    communityEvents.emplace_back(CommunityEvent::TIME_STEP);

    hasRun = true;
  }


  CKBDynamic::PowerlawCommunitySizeDistribution::PowerlawCommunitySizeDistribution(count minSize, count maxSize, double gamma, double alpha, double densityGamma, double minSplitRatio) : CommunitySizeDistribution(minSize, maxSize), sequence(minSize, maxSize, gamma), alpha(alpha), densityGamma(densityGamma), minSplitRatio(minSplitRatio) {
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
		//assert(getProbability(sizeA) == probabilityA);
		//assert(getProbability(sizeB) == probabilityB);

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
