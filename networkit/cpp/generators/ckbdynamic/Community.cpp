
#include "Community.h"
#include "CKBDynamicImpl.h"

namespace NetworKit {
	namespace CKBDynamicImpl {

		namespace {
			/**
			 * Returns number of steps you need to wait until the next success (edge) occurs.
			 */
			count get_next_edge_distance(const double log_cp, CKBDynamicImpl& generator) {
				return static_cast<count>(1 + floor(log(1.0 - generator.drawProbability()) / log_cp));
			}

		}

		Community::Community(const Community& o) : tlx::ReferenceCounter(), id(o.generator.nextCommunityId()), desiredSize(o.desiredSize), edges(o.edges), nonEdges(o.nonEdges), nodes(o.nodes), neighbors(o.neighbors), edgeProbability(o.edgeProbability), storeNonEdges(o.storeNonEdges), generator(o.generator), currentEvent(o.currentEvent) {

			for (auto e : edges) {
				generator.addEdge(e.first, e.second, false);
			}

			for (auto u : nodes) {
				generator.addNodeToCommunity(u, CommunityPtr(this));
				if (currentEvent != nullptr) {
					currentEvent->notifyNodeAddedToCommunity(u, CommunityPtr(this));
				}
			}

			generator.addCommunity(CommunityPtr(this));
		}


		Community::Community(CKBDynamicImpl& generator) : id(generator.nextCommunityId()), desiredSize(0), edgeProbability(0), storeNonEdges(edgeProbability > 0.6), generator(generator), currentEvent(nullptr) {
			//neighbors.min_load_factor(0.05);
			generator.addCommunity(CommunityPtr(this));
		}

		void Community::setDesiredNumberOfNodes(count size) {
			desiredSize = size;
			changeEdgeProbability(generator.communitySizeSampler->getCommunityDensity(size));
		}

		void Community::verifyInvariants() const {
			#ifndef NDEBUG
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
			#endif
		}

		void Community::removeEdge(node u, node v, bool nodeLeft, bool noEraseUNeighbors) {
			edges.erase(canonicalEdge(u, v));
			if (storeNonEdges) {
				nonEdges.insert(canonicalEdge(u, v));
				verifyInvariants();
			}
			neighbors[v].erase(u);
			if (!noEraseUNeighbors) {
				neighbors[u].erase(v);
			}

			generator.removeEdge(u, v, nodeLeft);
		}

		void Community::addEdge(node u, node v, bool nodeJoined) {
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

			generator.addEdge(u, v, nodeJoined);
		}

		std::pair<node, node> Community::edgeFromIndex(index i) const {
			node u = 1 + std::floor(-0.5 + std::sqrt(0.25 + 2.0 * i));
			node v = i - (u * (u - 1) / 2);

			return canonicalEdge(nodes.at(u), nodes.at(v));
		}


		void Community::removeNode(node u) {
			assert(nodes.contains(u)); // assert for gdb in gtest which catches exceptions
			if (!nodes.contains(u)) throw std::runtime_error("Node not in community!");

			for (node v : neighbors[u]) {
				removeEdge(u, v, true, true);
			}

			neighbors.erase(u);
			nodes.erase(u);

			if (currentEvent != nullptr) {
				currentEvent->notifyNodeRemovedFromCommunity(u, CommunityPtr(this));
			}
			generator.removeNodeFromCommunity(u, CommunityPtr(this));

			if (storeNonEdges) {
				for (node i = 0; i < nodes.size(); ++i) {
					nonEdges.erase(canonicalEdge(u, nodes.at(i)));
				}
				verifyInvariants();
			}
		}

		node Community::removeRandomNode() {
			assert(nodes.size() > 0); // assert for gdb in gtest which catches exceptions
			if (nodes.size() == 0) throw std::runtime_error("Error, no nodes in community!");

			const node u = nodes.at(generator.drawIndex(nodes.size()));
			removeNode(u);
			return u;
		}

		void Community::addNode(node u) {
			assert(!nodes.contains(u)); // assert for gdb in gtest which catches exceptions
			if (nodes.contains(u)) throw std::runtime_error("Node already in community!");
			// The global community has desired size 0
			assert(desiredSize == 0 || nodes.size() < desiredSize);
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

			neighbors.insert({u, robin_hood::unordered_flat_set<node>()});
			//auto nu = neighbors.try_emplace(u).first;
			//nu.value().min_load_factor(0.05);
			const double log_cp = std::log(1.0 - edgeProbability);

			for (node next = get_next_edge_distance(log_cp, generator) - 1; next < nodes.size(); next += get_next_edge_distance(log_cp, generator)) {
				const node v = nodes.at(next);
				if (u != v) {
					addEdge(u, nodes.at(next), true);
				}
			}

			if (currentEvent != nullptr) {
				currentEvent->notifyNodeAddedToCommunity(u, CommunityPtr(this));
			}
			generator.addNodeToCommunity(u, CommunityPtr(this));
		}

		count Community::drawDesiredNumberOfEdges(double prob) const {
			return generator.drawBinomial(getMaximumNumberOfEdges(), prob);
		}

		void Community::perturbEdges(double prob) {
			assert(prob <= 1);
			assert(prob >= 0);

			// Edges have probability 1 - nothing to perturb
			if (edgeProbability == 1.0) return;

			const count edgesToPerturb = generator.drawBinomial(edges.size(), prob);
			if (edgesToPerturb == 0) return;

			const count desiredNumberOfEdges = drawDesiredNumberOfEdges(edgeProbability);
			const count numEdges = edges.size();
			const count numNonEdges = getMaximumNumberOfEdges() - numEdges;

			// Community is clique and shall remain clique - nothing to perturb
			if (desiredNumberOfEdges == numEdges && numNonEdges == 0) return;

			count numEdgesToAdd = 0, numEdgesToRemove = 0;

			if (desiredNumberOfEdges >= numEdges) {
				numEdgesToAdd = edgesToPerturb;
				const count additionalEdgesDesired = desiredNumberOfEdges - numEdges;
				if (additionalEdgesDesired < edgesToPerturb) {
					numEdgesToRemove = edgesToPerturb - additionalEdgesDesired;
				}

				assert(edges.size() + numEdgesToAdd - numEdgesToRemove <= desiredNumberOfEdges);
			} else {
				assert(edges.size() > desiredNumberOfEdges);
				numEdgesToRemove = edgesToPerturb;
				const count edgesToLoose = numEdges - desiredNumberOfEdges;
				if (edgesToLoose < edgesToPerturb) {
					numEdgesToAdd = edgesToPerturb - edgesToLoose;
				}

				assert(edges.size() + numEdgesToAdd - numEdgesToRemove >= desiredNumberOfEdges);
			}

			assert(std::max(numEdgesToRemove, numEdgesToAdd) == edgesToPerturb);

			if (numEdgesToAdd > numNonEdges) {
				// As the desired number of edges can never be larger than
				// getMaximumNumberOfEdges() the following condition holds.
				assert(numEdgesToAdd - numNonEdges > numEdgesToRemove);

				numEdgesToRemove -= (numEdgesToAdd - numNonEdges);
				numEdgesToAdd = numNonEdges;
			}

			edges.random_sample(numEdgesToRemove);

			assert(edges.size() + numEdgesToAdd == getMaximumNumberOfEdges() || std::max(numEdgesToAdd, numEdgesToRemove) == edgesToPerturb);

			addRandomEdges(numEdgesToAdd);

			// remove non-edges if we will remove a lot of edges
			if (storeNonEdges && edges.size() - numEdgesToRemove < 0.25 * getMaximumNumberOfEdges()) {
				nonEdges.clear();
				storeNonEdges = false;
			}

			for (size_t i = numEdgesToRemove; i > 0; --i) {
				auto e = edges.at(i - 1);
				removeEdge(e.first, e.second, false);
			}

			verifyInvariants();
		}

		void Community::removeRandomEdges(count k) {
			assert(k <= edges.size());

			if (storeNonEdges && edges.size() - k < 0.25 * getMaximumNumberOfEdges()) {
				nonEdges.clear();
				storeNonEdges = false;
			}

			for (count i = 0; i < k; ++i) {
				auto e = edges.at(generator.drawIndex(edges.size()));
				removeEdge(e.first, e.second, false);
			}
		}

		void Community::addRandomEdges(count k) {
			const count numEdgesWanted = edges.size() + k;
			assert(numEdgesWanted <= getMaximumNumberOfEdges());

			if (!storeNonEdges && numEdgesWanted >= 0.75 * getMaximumNumberOfEdges()) {
				nonEdges.reserve(getMaximumNumberOfEdges() - edges.size());
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

			if (storeNonEdges) {
				assert(k <= nonEdges.size());
				verifyInvariants();

				for (count i = 0; i < k; ++i) {
					auto e = nonEdges.at(generator.drawIndex(nonEdges.size()));
					addEdge(e.first, e.second, false);
				}
			} else {
				const count n = nodes.size();

				while (edges.size() < numEdgesWanted) {
					node u = nodes.at(generator.drawIndex(n));
					node v = nodes.at(generator.drawIndex(1, n));
					if (u == v) v = nodes.at(0);
					assert(u != v);
					if (!edges.contains(canonicalEdge(u, v))) {
						addEdge(u, v, false);
					}
				}
			}
		}

		void Community::changeEdgeProbability(double prob) {
			if (prob > 1) prob = 1;

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

			while (other.edges.size() > 0) {
				const auto e = other.edges.at(0);
				// First add edge here to ensure it exists at least once globally
				// so we don't generate remove/add events globally.
				if (!edges.contains(e)) {
					addEdge(e.first, e.second, false);
				}

				other.removeEdge(e.first, e.second, false);
			}

			assert(other.edges.size() == 0);

			while (other.nodes.size() > 0) {
				const node u = other.nodes.at(0);
				other.removeNode(u);
			}

			assert(other.neighbors.size() == 0);
			assert(other.nonEdges.size() == 0);
			assert(other.nodes.size() == 0);
		}

		void Community::setCurrentEvent(CommunityChangeEvent *event) {
			if (event != this->currentEvent) {
				this->currentEvent = event;
				generator.addCommunity(CommunityPtr(this));
			}
		}

		bool Community::canRemoveNode() const {
			return currentEvent == nullptr || currentEvent->canRemoveNode();
		}
	}
}
