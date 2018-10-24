#ifndef CKBDYNAMIC_BUCKETSAMPLING_H_
#define CKBDYNAMIC_BUCKETSAMPLING_H_

#include "../../Globals.h"
#include "../PowerlawDegreeSequence.h"
#include "../../auxiliary/SamplingSet.h"

namespace NetworKit {
	namespace CKBDynamicImpl {
		class BucketSampling {
		private:
			static constexpr count oversample_fraction = 5;

			struct node_data {
				// positions in the full_slots vector where this node occurs
				std::vector<index> full_slot_positions;
				// fractional part of the number of slots times oversample_fraction
				// this is negative if the node is in too many communities
				int64_t current_fractional_slot;
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
			BucketSampling(count n, count minSlots, count maxSlots, double exponent);

			/**
			 * Sample @a communitySize nodes from the available memberships.
			 *
			 * @note This does not modify the distribution. Run
			 * BucketSampling::assignCommunity() for every node to actually take the slots.
			 * @param communitySize The number of nodes to sample.
			 * @return The sampled nodes.
			 */
			std::vector<node> birthCommunityNodes(count communitySize, const Aux::SamplingSet<node>& existingNodes = Aux::SamplingSet<node>());

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
			node removeNode(node nodeId, bool withFraction);
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
			node removeNode(node nodeId);

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
	}
}
#endif
