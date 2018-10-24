
#include "Community.h"
#include "CKBDynamicImpl.h"
#include "../../auxiliary/UniqueSampler.h"

namespace NetworKit {
	namespace CKBDynamicImpl {

	namespace {
		/**
		 * Returns number of steps you need to wait until the next success (edge) occurs.
		 */
		count get_next_edge_distance(const double log_cp) {
			return static_cast<count>(1 + floor(log(1.0 - Aux::Random::probability()) / log_cp));
		}

	}

	Community::Community(const Community& o) : id(o.generator.nextCommunityId()), edges(o.edges), nonEdges(o.nonEdges), nodes(o.nodes), neighbors(o.neighbors), edgeProbability(o.edgeProbability), available(o.available), storeNonEdges(o.storeNonEdges), generator(o.generator) {
		for (auto e : edges) {
			generator.addEdge(e.first, e.second);
		}

		for (auto u : nodes) {
			generator.addNodeToCommunity(u, CommunityPtr(this));
		}

		generator.addCommunity(CommunityPtr(this));
	}


	Community::Community(double edgeProbability, CKBDynamicImpl& generator) : id(generator.nextCommunityId()), edgeProbability(edgeProbability), available(false), storeNonEdges(edgeProbability > 0.6), generator(generator) {
		generator.addCommunity(CommunityPtr(this));
	}

	void Community::verifyInvariants() const {
		if (storeNonEdges) {
			for (auto e : edges) {
				assert(generator.canonicalEdge(e.first, e.second) == e);
				assert(!nonEdges.contains(e));
				assert(nodes.contains(e.first));
				assert(nodes.contains(e.second));
			}

			for (auto e : nonEdges) {
				assert(generator.canonicalEdge(e.first, e.second) == e);
				assert(nodes.contains(e.first));
				assert(nodes.contains(e.second));
			}

			assert(nonEdges.size() + edges.size() == getMaximumNumberOfEdges());
		}
	}

	void Community::removeEdge(node u, node v) {
			edges.erase(generator.canonicalEdge(u, v));
			if (storeNonEdges) {
	nonEdges.insert(generator.canonicalEdge(u, v));
				verifyInvariants();
			}
			neighbors[v].erase(u);
			neighbors[u].erase(v);
			generator.removeEdge(u, v);
	}

	void Community::addEdge(node u, node v) {
		assert(u != v);
		assert(nodes.contains(u));
		assert(nodes.contains(v));
			neighbors[u].insert(v);
			neighbors[v].insert(u);
			edges.insert(generator.canonicalEdge(u, v));
			if (storeNonEdges) {
	nonEdges.erase(generator.canonicalEdge(u, v));
				verifyInvariants();
			}
			generator.addEdge(u, v);
	}

	std::pair<node, node> Community::edgeFromIndex(index i) const {
		node u = 1 + std::floor(-0.5 + std::sqrt(0.25 + 2.0 * i));
		node v = i - (u * (u - 1) / 2);

		return generator.canonicalEdge(nodes.at(u), nodes.at(v));
	}


	void Community::removeNode(node u) {
		assert(nodes.contains(u));
		if (!nodes.contains(u)) throw std::runtime_error("Node not in community!");

		for (node v : neighbors[u]) {
			// We cannot delete from neighbors[u] because we iterate over it!
			edges.erase(generator.canonicalEdge(u, v));
			neighbors[v].erase(u);
			generator.removeEdge(u, v);
		}

		neighbors.erase(u);
		nodes.erase(u);
		generator.removeNodeFromCommunity(u, CommunityPtr(this));

		if (storeNonEdges) {
			for (node i = 0; i < nodes.size(); ++i) {
	nonEdges.erase(generator.canonicalEdge(u, nodes.at(i)));
			}
			verifyInvariants();
		}
	}

	void Community::removeRandomNode() {
		assert(nodes.size() > 0);
		if (nodes.size() == 0) throw std::runtime_error("Error, no nodes in community!");

		node u = nodes.at(Aux::Random::index(nodes.size()));
		removeNode(u);
	}

	void Community::addNode(node u) {
		assert(!nodes.contains(u));
		if (nodes.contains(u)) throw std::runtime_error("Node already in community!");
		nodes.insert(u);

		if (storeNonEdges) {
			for (node i = 0; i < nodes.size(); ++i) {
				const node v = nodes.at(i);
				if (u != v) {
					nonEdges.insert(generator.canonicalEdge(u, v));
				}
			}

			verifyInvariants();
		}

		neighbors.insert(std::make_pair(u, Aux::SamplingSet<node>()));
		const double log_cp = std::log(1.0 - edgeProbability);

		for (node next = get_next_edge_distance(log_cp) - 1; next < nodes.size(); next += get_next_edge_distance(log_cp)) {
			const node v = nodes.at(next);
			if (u != v) {
				addEdge(u, nodes.at(next));
			}
		}

		generator.addNodeToCommunity(u, CommunityPtr(this));
	}

	count Community::drawDesiredNumberOfEdges(double prob) const {
			std::binomial_distribution<count> distr(getMaximumNumberOfEdges(), prob);
			return distr(Aux::Random::getURNG());
	}

	void Community::perturbEdges(double prob) {
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

	void Community::removeRandomEdges(count k) {
		assert(k <= edges.size());

		for (count i = 0; i < k; ++i) {
			auto e = edges.at(Aux::Random::index(edges.size()));
			removeEdge(e.first, e.second);
		}
	}

	void Community::addRandomEdges(count k) {
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
	if (u != v && !edges.contains(generator.canonicalEdge(u, v))) {
		addEdge(u, v);
	}
			}
		}
	}

	void Community::changeEdgeProbability(double prob) {
		if (prob > 1) prob = 1;

		if (prob < 0.4 && storeNonEdges) {
			nonEdges.clear();
			storeNonEdges = false;
		} else if (prob > 0.6 && !storeNonEdges) {

			for (index i = 0; i < nodes.size(); ++i) {
	const node u = nodes.at(i);
	for (index j = i + 1; j < nodes.size(); ++j) {
		const node v = nodes.at(j);

		const auto e = generator.canonicalEdge(u, v);
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

	void Community::combineWith(Community& other) {
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

	void Community::setAvailable( bool avail ) {
		if (avail) {
			assert(nodes.size() > 0);
		}
		if (this->available != avail) {
			this->available = avail;
			generator.addCommunity(CommunityPtr(this));
		}
	}
	}
}
