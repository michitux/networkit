#ifndef CKB_DYNAMIC_H
#define CKB_DYNAMIC_H

#include "../graph/Graph.h"
#include "../structures/Cover.h"
#include "../dynamics/GraphEvent.h"
#include "../dynamics/CommunityEvent.h"
#include "../base/Algorithm.h"
#include "PowerlawDegreeSequence.h"
#include <tlx/container/btree_set.hpp>
#include <tlx/container/btree_map.hpp>
#include <tlx/counting_ptr.hpp>
#include "../auxiliary/SamplingSet.h"

namespace NetworKit {
	class CKBDynamic : public Algorithm {
		CKBDynamic(const Cover& model);
		CKBDynamic(count n, const Cover& model);
		CKBDynamic(count n, count minCommunitySize, count maxCommunitySize, count minCommunityMembership, count maxCommunityMembership, double eventProbability, double intraCommunityEdgeProbability, double epsilon);

		virtual void run() override;

		std::vector<GraphEvent> getGraphEvents();
		std::vector<CommunityEvent> getCommunityEvents();
	private:
		std::vector<GraphEvent> graphEvents;
		std::vector<CommunityEvent> communityEvents;

		class Community : public tlx::ReferenceCounter {
		public:
			/**
			 * Create an exact copy of a given community @a o
			 *
			 * @param o The community to copy.
			 */
			Community(const Community& o);

			Community(double edgeProbability, CKBDynamic& generator);

			/**
			 * Remove the given node @a u from the community.
			 *
			 * This removes all edges incident to @a u from the community.
			 */
			void removeNode(node u);

			/**
			 * Adds a node to the community and generates
			 * edges from @a u to the community with the
			 * current edge probability of the community.
			 *
			 * @param u The node to add.
			 */
			void addNode(node u);

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

			count getNumberOfEdges() const { return edges.size(); };
		private:
			void removeEdge(node u, node v);
			void addEdge(node u, node v);

			Aux::SamplingSet<std::pair<node, node>> edges;
			// only used if edgeProbability > 0.5.
			Aux::SamplingSet<std::pair<node, node>> nonEdges;
			Aux::SamplingSet<node> nodes;
			tlx::btree_map<node, Aux::SamplingSet<node>> neighbors;
			double edgeProbability;
			bool isAvailable;
			bool storeNonEdges;
			CKBDynamic& generator;
		};

		class CommunityChangeEvent {
		public:
			CommunityChangeEvent(CKBDynamic& generator);
			virtual void nextStep() = 0;
			bool isActive();
		private:
			bool active;
			count numSteps;
			count currentStep;
			CKBDynamic& generator;
		};

		// Create community from a subset of given nodes.
		// Adds further nodes over time.
		class CommunityBirthEvent : public CommunityChangeEvent {
		public:
			CommunityBirthEvent(Community &community, std::vector<node> nodes, count numSteps);
		private:
			std::vector<node> nodes;
			Community& community;
		};

		// Remove nodes over time.
		// Lets community die at the end.
		class CommunityDeathEvent : public CommunityChangeEvent {
		public:
			CommunityDeathEvent(Community &community, count numSteps);
		private:
			Community& community;
		};

		// Takes two communities.
		// Add nodes from community A to B
		// and vice-versa. However, while adding nodes,
		// decrease the edge probability until it is so low
		// that when combining both communities, we get the
		// target edge probability (note: an existing overlap
		// needs to be considered here). In the final step,
		// just combine both edge sets and remove duplicate
		// edges (needs special support in the community
		// object).
		class CommunityMergeEvent : public CommunityChangeEvent {
		public:
			CommunityMergeEvent(Community& communityA, Community& communityB, double targetEdgeProbability);
		private:
			Community &communityA;
			Community &communityB;
		};

		// Takes one community.
		// Creates a second community as exact copy with the same edge set.
		// Removes parts of both communities over time.
		// Increases the density of both communities over time.
		class CommunitySplitEvent : public CommunityChangeEvent {
		};


		class CommunitySizeDistribution {
		public:
			virtual std::pair<count, double> drawCommunity() = 0;
			virtual std::pair<count, double> mergeCommunities(count sizeA, double probabilityA, count sizeB, double probabilityB, count combinedNodes) = 0;
			virtual std::pair<std::pair<count, double>, std::pair<count, double>> splitCommunity(count size, double probability) = 0;
		};

		class PowerlawCommunitySizeDistribution : public CommunitySizeDistribution {
		private:
			PowerlawDegreeSequence sequence;
			double alpha;
			double densityGamma;
			double minSplitRatio;

			double getProbability(count size) const;
		public:
			PowerlawCommunitySizeDistribution(count minSize, count maxSize, double gamma, double alpha, double densityGamma, double minSplitRatio);

			virtual std::pair<count, double> drawCommunity() override;

			virtual std::pair<count, double> mergeCommunities(count sizeA, double probabilityA, count sizeB, double probabilityB, count combinedNodes) override;

			virtual std::pair<std::pair<count, double>, std::pair<count, double>> splitCommunity(count size, double probability) override;
		};

		class BucketSampling {
		private:
			static constexpr count oversample_fraction = 5;

			struct node_data {
				// positions in the full_slots vector where this node occurs
				std::vector<index> full_slot_positions;
				// fractional part of the number of slots times oversample_fraction
				// this is negative if the node is in too many communities
				index current_fractional_slot;
				// position in the fractional slot vector of current_fractional_slot > 0
				index fractional_slot_position;
				// desired degree of the node. 0 if the node has been removed.
				count degree;
				// position in the nodes_by_memberships vector if degree > 0.
				index nodes_by_memberships_pos;
			};

			struct slot {
				node u; /* The node for which the slot is */
				index pos_at_node; /* The position in the node's slot position vector that links to the slot */
			};

			PowerlawDegreeSequence pl_generator;
			// all nodes sorted by the number of desired memberships
			std::vector<std::vector<node>> nodes_by_memberships;
			std::vector<std::array<node, oversample_fraction-1>> nodes_with_fractional_slots_for_degree;
			std::vector<node_data> nodes;
			std::vector<bool> node_sampled_in_current_call;
			std::vector<slot> full_slots;
			std::array<std::vector<node>, oversample_fraction> fractional_slots;
			std::mt19937_64 random_engine;

			count sumOfDesiredMemberships;

			void verifyInvariants(bool verify_oversampling = false) const;
		public:
			/**
			 * Initialize the sampling data structure.
			 *
			 * @param n The number of initial nodes.
			 * @param minSlots The minimum number of memberships per node.
			 * @param maxSlots The maximum number of memberships per node.
			 * @param exponent The powerlaw exponent of the membership distribution.
			 * @param seed Seed for the random number generator used for various random decisions.
			 */
			BucketSampling(count n, count minSlots, count maxSlots, double exponent, uint64_t seed);

			/**
			 * Sample @a communitySize nodes from the available memberships.
			 *
			 * @note This does not modify the distribution. Run
			 * BucketSampling::assignCommunity() for every node to actually take the slots.
			 * @param communitySize The number of nodes to sample.
			 * @return The sampled nodes.
			 */
			std::vector<node> birthCommunityNodes(count communitySize);

			/**
			 * Assign the node @a nodeId to a community, i.e., take a slot from it.
			 *
			 * @param nodeId The node to assign.
			 */
			void assignCommunity(node nodeId);

			/**
			 * Remove the node @a nodeId from a community, i.e., assign a slot to it.
			 *
			 * This has no effect if the node has already been deleted, i.e.,
			 * if its desired number of memberships is 0.
			 *
			 * @param nodeId The node to remove from a community.
			 */
			void leaveCommunity(node nodeId);

		private:
			void removeFraction(node nodeId, count fraction);

			void addFraction(node nodeId, count fraction);

			node addNode(count degree, bool oversample);

		public:

			/**
			 * Add a new node with a randomly sampled number of desired memberships.
			 *
			 * This adds exactly 1.2 times the number of desired
			 * memberships slots to the sampling data structure. For every
			 * 5 nodes of the same number of desired membership, an
			 * additional 6th node is added internally, so this may
			 * actually add two nodes.
			 *
			 * @return The id of the added node.
			 */
			node addNode();

			/**
			 * Add a new node with the given number of desired memberships.
			 *
			 * This adds exactly 1.2 times the number of desired
			 * memberships slots to the sampling data structure. For every
			 * 5 nodes of the same number of desired membership, an
			 * additional 6th node is added internally, so this may
			 * actually add two nodes.
			 *
			 * @param degree The number of desired memberships of the node to add.
			 * @return The id of the added node.
			 */
			node addNode(count degree);

		private:
			void removeNode(node nodeId, bool withFraction);
		public:
			/**
			 * Remove a node from the sampling data structure.
			 *
			 * This removes exactly 1.2 times the number of desired
			 * memberships of the node slots from the sampling data
			 * structure. This may remove an additional node with the same
			 * number of desired memberships. This does nothing if the
			 * node has already been deleted.
			 *
			 * @param nodeId The id of the node to remove.
			 */
			void removeNode(node nodeId);

			/**
			 * Get the total number of nodes. Some of these nodes maybe deleted.
			 *
			 * All node ids are between 0 and getNumberOfNodes()-1.
			 *
			 * @return The total number of nodes.
			 */
			count getNumberOfNodes() const {
				return nodes.size();
			}

			/**
			 * Get the sum of desired memberships.
			 *
			 * This only includes explicitly added nodes and not nodes
			 * added due to the oversampling.
			 *
			 * @return The sum of desired memberships.
			 */
			count getSumOfDesiredMemberships() const {
				return sumOfDesiredMemberships;
			}

			/**
			 * Get the number of desired memberships of a node.
			 *
			 * This does not include oversampling. However, if due to the
			 * oversampling additional nodes are added, they will also
			 * report a positive number of desired memberships even though
			 * they are not counted for @a getSumOfDesiredMemberships().
			 * This reports 0 for deleted nodes.
			 *
			 * @param nodeId The id of the node.
			 * @return The number of desired memberships of the node.
			 */
			count getDesiredMemberships(count nodeId) const {
				return nodes[nodeId].degree;
			}

			void printSlots() const {
				std::vector<count> slotsPerDegree(pl_generator.getMaximumDegree() + 1);

				for (const node_data& node : nodes) {
					slotsPerDegree[node.degree] += node.full_slot_positions.size();
				}

				for (count d = 0; d < pl_generator.getMaximumDegree() + 1; ++d) {
					if (slotsPerDegree[d] > 0) {
						std::cout << "degree " << d << " has " << slotsPerDegree[d] << " slots" << std::endl;
					}
				}
			}
		};

		using CommunityPtr = tlx::CountingPtr<Community>;

		void addEdge(node u, node v);
		void removeEdge(node u, node v);
		void addNodeToCommunity(node u, CommunityPtr com);
		void removeNodeFromCommunity(node u, CommunityPtr com);
		void addAvailableCommunity(CommunityPtr com);

		Aux::SamplingSet<CommunityPtr> splittableCommunities;
		Aux::SamplingSet<CommunityPtr> availableCommunities;
		std::vector<tlx::btree_set<CommunityPtr>> nodeCommunities;
	};

} // namespace NetworKit

#endif
