#ifndef CKBDYNAMIC_COMMUNITY_H_
#define CKBDYNAMIC_COMMUNITY_H_

#include <tlx/counting_ptr.hpp>
#include <tlx/math/integer_log2.hpp>
#include "../../Globals.h"
#include "../../auxiliary/SamplingSet.h"
#include "NodePairHash.h"
#include <tsl/robin_map.h>

namespace NetworKit {
	namespace CKBDynamicImpl {
		class CKBDynamicImpl;
		class CommunityChangeEvent;

		class Community : public tlx::ReferenceCounter {
		public:
			static std::pair<node, node> canonicalEdge(node u, node v) {
				if (u < v) return std::make_pair(u, v);
				return std::make_pair(v, u);
			}
			/**
			 * Create an exact copy of a given community @a o
			 *
			 * @param o The community to copy.
			 */
			Community(const Community& o);

			Community(CKBDynamicImpl& generator);

			/**
			 * Remove the given node @a u from the community.
			 *
			 * This removes all edges incident to @a u from the community.
			 */
			void removeNode(node u);

			/**
			 * Remove a random node from the community.
			 *
			 * @return The removed node.
			 */
			node removeRandomNode();

			/**
			 * Adds a node to the community and generates
			 * edges from @a u to the community with the
			 * current edge probability of the community.
			 *
			 * @param u The node to add.
			 */
			void addNode(node u);

			double getEdgeProbability() const { return edgeProbability; }

			/**
			 * Change the probability of edges to @a prob.
			 *
			 * This adds or removes random edges to achieve the desired edge probability.
			 *
			 * @param u The node to remove.
			 */
			void changeEdgeProbability(double prob);

			/**
			 * Combine the edges of this community with the edges of another community.
			 * This assumes that the nodes sets are identical.
			 * Overlapping edges are removed and the other community will be empty afterwards.
			 *
			 * @param other The community to combine the nodes with.
			 */
			void combineWith(Community& other);

			/**
			 * Perturb edges with a certain probability. 0
			 * changes no edges, 1 removes all edges and
			 * completely re-generates the community.
			 *
			 * FIXME possibly we want to call this with an number of edges to perturb.
			 *
			 * @param p The probability/percentage with which edges shall be removed or added.
			 */
			void perturbEdges(double prob);


			count getNumberOfNodes() const { return nodes.size(); };

			count getDesiredNumberOfNodes() const { return desiredSize; };

			/**
			 * Set the desired number of nodes.
			 *
			 * This automatically also adjusts the edge probability using the community size sampler.
			 *
			 * @param size The desired size
			 */
			void setDesiredNumberOfNodes(count size);

			count getNumberOfEdges() const { return edges.size(); };

			count getMaximumNumberOfEdges() const { return (nodes.size() * (nodes.size() - 1) / 2); }

			const Aux::SamplingSet<node>& getNodes() const { return nodes; };

			bool hasNode(node u) const { return nodes.contains(u); };

			index getId() const { return id; };

			bool isAvailable() const { return currentEvent == nullptr; }

			void setCurrentEvent(CommunityChangeEvent* event);

			bool canRemoveNode() const;
		private:
			void removeEdge(node u, node v, bool nodeLeft);
			void addEdge(node u, node v, bool nodeJoined);

			void removeRandomEdges(count k);
			void addRandomEdges(count k);

			count drawDesiredNumberOfEdges(double prob) const;

			void verifyInvariants() const;

			std::pair<node, node> edgeFromIndex(index i) const;

			index id;
			count desiredSize;
			Aux::SamplingSet<std::pair<node, node>, NodePairHash> edges;
			// only used if edgeProbability > 0.5.
			Aux::SamplingSet<std::pair<node, node>, NodePairHash> nonEdges;
			Aux::SamplingSet<node> nodes;
			tsl::robin_map<node, Aux::SamplingSet<node>> neighbors;
			double edgeProbability;
			bool storeNonEdges;
			CKBDynamicImpl& generator;
			CommunityChangeEvent* currentEvent;
		};

		using CommunityPtr = tlx::CountingPtr<Community>;
	}
}

namespace std {
	// inspired by https://stackoverflow.com/questions/20953390/what-is-the-fastest-hash-function-for-pointers
	template <>
	struct hash<NetworKit::CKBDynamicImpl::CommunityPtr> {
		std::size_t operator()(const NetworKit::CKBDynamicImpl::CommunityPtr& k) const {
			using std::size_t;
			static const size_t shift = tlx::integer_log2_floor(1 + sizeof(NetworKit::CKBDynamicImpl::Community));
			return reinterpret_cast<size_t>(k.get()) >> shift;
		};
	};
}

#endif
