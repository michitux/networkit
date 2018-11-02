#include "BucketSampling.h"
#include <tlx/unused.hpp>
#include "../../auxiliary/UniqueSampler.h"

namespace NetworKit {
	namespace CKBDynamicImpl {
		void BucketSampling::verifyInvariants(bool verify_oversampling) const {
			tlx::unused(verify_oversampling);
#ifndef NDEBUG

			count num_full_slots = 0;
			std::array<count, oversample_fraction> num_fractional_slots;
			num_fractional_slots.fill(0);

			count num_existing_nodes = 0;

			for (count u = 0; u < nodes.size(); ++u) {
				const node_data& node = nodes[u];

				if (node.degree == 0) {
					assert(node.full_slot_positions.empty());
					assert(node.current_fractional_slot <= 0);
				} else {
					++num_existing_nodes;

					if (node.current_fractional_slot < 0) {
						assert(node.full_slot_positions.empty());
					}

					for (index i = 0; i < node.full_slot_positions.size(); ++i) {
						const index slot_pos = node.full_slot_positions[i];
						++num_full_slots;

						assert(full_slots[slot_pos].u == u);
						assert(full_slots[slot_pos].pos_at_node == i);
					}

					if (node.current_fractional_slot > 0) {
						++num_fractional_slots[node.current_fractional_slot];
						assert(fractional_slots[node.current_fractional_slot][node.fractional_slot_position] == u);

						auto& nodes_with_frac = nodes_with_fractional_slots_for_degree[node.degree];
						assert(std::find(nodes_with_frac.begin(), nodes_with_frac.end(), u) != nodes_with_frac.end());
					}

					assert(nodes_by_memberships[node.degree][node.nodes_by_memberships_pos] == u);
				}
			}

			assert(num_full_slots == full_slots.size());
			for (count i = 0; i < oversample_fraction; ++i) {
				assert(num_fractional_slots[i] == fractional_slots[i].size());
			}

			count membership_sum = 0;
			for (count i = 0; i <= pl_generator.getMaximumDegree(); ++i) {
				if (i < pl_generator.getMinimumDegree()) {
					assert(nodes_by_memberships[i].empty());
					for (node v : nodes_with_fractional_slots_for_degree[i]) {
						assert(v == none);
					}
				} else if (verify_oversampling) {
					node nodes_found = 0;

					for (node v : nodes_with_fractional_slots_for_degree[i]) {
						if (v != none) {
							assert(nodes[v].degree == i);
							++nodes_found;
						}
					}

					assert(nodes_found == (nodes_by_memberships[i].size() % (oversample_fraction + 1)));
				}

				membership_sum += nodes_by_memberships[i].size();
			}

			assert(membership_sum == num_existing_nodes);
#endif
		}


		BucketSampling::BucketSampling(count n, count minSlots, count maxSlots, double exponent) : pl_generator(minSlots, maxSlots, exponent), nodes_by_memberships(maxSlots + 1), nodes_with_fractional_slots_for_degree(maxSlots + 1), sumOfDesiredMemberships(0) {
			pl_generator.run();

			for (count d = 0; d <= maxSlots; ++d) {
				nodes_with_fractional_slots_for_degree[d].fill(none);
			}

			for (count u = 0; u < n; ++u) {
				addNode();
			}
		}

		std::vector<node> BucketSampling::birthCommunityNodes(count communitySize, const Aux::SamplingSet<node>& existingNodes) {
			std::vector<node> result;

			std::array<count, oversample_fraction> slot_boundaries;

			slot_boundaries[0] = full_slots.size() * oversample_fraction;
			for (count i = 1; i < oversample_fraction; ++i) {
				// The weight of each slot in this array
				const count total_slot_weight = fractional_slots[i].size() * i;
				slot_boundaries[i] = slot_boundaries[i-1] + total_slot_weight;
			}

			const count max_slot = slot_boundaries.back();
			Aux::UniqueSampler slotSampler(max_slot);

			// Give up after trying every slot once.
			for (count attempt = 0; attempt < max_slot && result.size() < communitySize; ++attempt) {
				const count sampled_slot = slotSampler.draw();
				node u = 0;

				if (sampled_slot < slot_boundaries[0]) {
					u = full_slots[sampled_slot / oversample_fraction].u;
				} else {
					for (count i = 1; i < oversample_fraction; ++i) {
						if (sampled_slot < slot_boundaries[i]) {
							const count slot_in_bin = sampled_slot - slot_boundaries[i - 1];
							const count index_in_bin = slot_in_bin / i;
							u = fractional_slots[i][index_in_bin];
							break;
						}
					}
				}

				if (!node_sampled_in_current_call[u] && !existingNodes.contains(u)) {
					result.push_back(u);
					node_sampled_in_current_call[u] = true;
				}
			}

			for (node u : result) {
				node_sampled_in_current_call[u] = false;
			}

			return result;
		}

		void BucketSampling::assignCommunity(node nodeId) {
			if (!nodes[nodeId].full_slot_positions.empty()) {
				// relocate slot to the end of full_slots
				index slot_pos = nodes[nodeId].full_slot_positions.back();
				assert(full_slots[slot_pos].u == nodeId);

				std::swap(full_slots[slot_pos], full_slots.back());

				// repair the reference to the other slot position
				slot& other_slot = full_slots[slot_pos];
				assert(nodes[other_slot.u].full_slot_positions[other_slot.pos_at_node] == full_slots.size() - 1);
				nodes[other_slot.u].full_slot_positions[other_slot.pos_at_node] = slot_pos;
				nodes[nodeId].full_slot_positions.pop_back();
				full_slots.pop_back();
				verifyInvariants();
			} else {
				removeFraction(nodeId, oversample_fraction);
				verifyInvariants();
			}

		}

		void BucketSampling::leaveCommunity(node nodeId) {
			node_data& node = nodes[nodeId];

			/* The node was removed internally but was still in a community. Do not give it back any slots. */
			if (node.degree == 0) return;

			if (node.current_fractional_slot < 0) {
				assert(node.full_slot_positions.empty());

				addFraction(nodeId, oversample_fraction);
			} else {
				index slot_pos = node.full_slot_positions.size();
				node.full_slot_positions.push_back(full_slots.size());
				full_slots.push_back({nodeId, slot_pos});
			}

			verifyInvariants();
		}

		void BucketSampling::removeFraction(node nodeId, count fraction) {
			if (fraction == 0) return;

			node_data& u = nodes[nodeId];

			if (u.current_fractional_slot > 0) {
				std::vector<count>& frac_slots = fractional_slots[u.current_fractional_slot];
				// relocate slot to the end of the fractional slots
				index slot_pos = u.fractional_slot_position;
				std::swap(frac_slots[slot_pos], frac_slots.back());

				node other_node = frac_slots[slot_pos];
				// repair the reference to the other slot position
				nodes[other_node].fractional_slot_position = slot_pos;
				// remove the slot
				frac_slots.pop_back();
			}

			u.current_fractional_slot -= fraction;
		}

		void BucketSampling::addFraction(node nodeId, count fraction) {
			if (fraction == 0) return;

			node_data& u = nodes[nodeId];

			u.current_fractional_slot += fraction;

			if (u.current_fractional_slot > 0) {
				std::vector<node>& frac_slots = fractional_slots[u.current_fractional_slot];
				u.fractional_slot_position = frac_slots.size();
				frac_slots.push_back(nodeId);
			}
		}

		node BucketSampling::addNode(count degree, bool oversample) {
			const node u = nodes.size();
			nodes.emplace_back();
			node_sampled_in_current_call.push_back(false);

			nodes[u].degree = degree;
			nodes[u].nodes_by_memberships_pos = nodes_by_memberships[degree].size();
			nodes_by_memberships[degree].push_back(u);

			int64_t u_fractional_slot = degree % oversample_fraction;
			count u_additional_full_slots = degree / oversample_fraction;
			count num_full_slots = degree;

			if (oversample) {
				auto& nodes_with_fractional_slots = nodes_with_fractional_slots_for_degree[degree];
				if (nodes_by_memberships[degree].size() % (oversample_fraction + 1) == oversample_fraction) {
					// remove oversampled slots
					for (node& v : nodes_with_fractional_slots) {
						assert(v != none);

						for (count i = 0; i < u_additional_full_slots; ++i) {
							assignCommunity(v);
						}

						removeFraction(v, u_fractional_slot);
						v = none;
					}

					// create a new node with exactly degree slots
					addNode(degree, false);

					// no fractional slot for this node
					u_fractional_slot = 0;
				} else {
					num_full_slots += u_additional_full_slots;
					addFraction(u, u_fractional_slot);

					bool found = false;
					for (node& v : nodes_with_fractional_slots) {
						if (v == none) {
							v = u;
							found = true;
							break;
						}
					}

					assert(found);
				}
			} else {
				assert(nodes_by_memberships[degree].size() % (oversample_fraction + 1) == 0);
				u_fractional_slot = 0;
			}

			nodes[u].full_slot_positions.reserve(num_full_slots);

			for (count i = 0; i < num_full_slots; ++i) {
				leaveCommunity(u);
			}

			assert(nodes[u].full_slot_positions.size() == num_full_slots);
			assert(nodes[u].current_fractional_slot == u_fractional_slot);

			return u;
		}

		node BucketSampling::addNode() {
			return addNode(pl_generator.getDegree());
		}

		node BucketSampling::addNode(node degree) {
			node u = addNode(degree, true);
			sumOfDesiredMemberships += degree;
			verifyInvariants(true);
			return u;
		}

		node BucketSampling::removeNode(node nodeId, bool withFraction) {
			node additionallyRemovedNode = none;
			// make sure there are no more slots for this node
			while (nodes[nodeId].current_fractional_slot >= 0) {
				assignCommunity(nodeId);
			}

			const count degree = nodes[nodeId].degree;
			nodes[nodeId].degree = 0;

			{ // delete node from nodes_by_memberships,
				const node partner = nodes_by_memberships[degree].back();
				std::swap(
					  nodes_by_memberships[degree][nodes[nodeId].nodes_by_memberships_pos],
					  nodes_by_memberships[degree].back()
					  );
				nodes[partner].nodes_by_memberships_pos = nodes[nodeId].nodes_by_memberships_pos;
				nodes_by_memberships[degree].pop_back();
			}

			auto& nodes_with_fractional_slots = nodes_with_fractional_slots_for_degree[degree];
			bool hadFraction = false;
			count nodesWithFraction = 0;
			for (node& v : nodes_with_fractional_slots) {
				if (v == nodeId) {
					v = none;
					hadFraction = true;
				} else if (v != none) {
					++nodesWithFraction;
				}
			}

			assert(withFraction || !hadFraction);

			// if node had a fraction, we are
			// Otherwise, we need to get the fractional part from somewhere else.
			if (!hadFraction && withFraction) {
				count u_fractional_slot = degree % oversample_fraction;
				count u_additional_full_slots = degree / oversample_fraction;

				if (nodesWithFraction > 0) {
					// if there is another node with a fractional part, it is easy, too:
					// just remove the additional slots from this node.
					std::uniform_int_distribution<index> random_selector(0, nodesWithFraction-1);
					index pos = random_selector(Aux::Random::getURNG());

					index i = 0;
					for (node& v : nodes_with_fractional_slots) {
						if (v != none) {
							if (i == pos) {
								// remove additional fractional part from this node
								for (index j = 0; j < u_additional_full_slots; ++j) {
									assignCommunity(j);
								}

								removeFraction(v, u_fractional_slot);

								v = none;
								break;
							}
							++i;
						}
					}
				} else {
					// now the difficult case: we need to
					// a) remove a complete other node and
					// b) find oversample_fraction-1 nodes of the same degree that then get a fraction!
					assert(nodes_by_memberships[degree].size() >= oversample_fraction);
					std::vector<node> sampled_nodes;

					// if we have exactly as many nodes as we want, just shuffle them
					if (nodes_by_memberships[degree].size() == oversample_fraction) {
						sampled_nodes = nodes_by_memberships[degree];
						std::shuffle(sampled_nodes.begin(), sampled_nodes.end(), Aux::Random::getURNG());
					} else {
						// now we should have at least twice as many nodes
						// as we need, so in expectation at least every second
						// sampled node has not been sampled yet.
						assert(nodes_by_memberships[degree].size() >= 2 * oversample_fraction);
						sampled_nodes.reserve(oversample_fraction);
						std::uniform_int_distribution<int> sample(0, nodes_by_memberships[degree].size() - 1);
						while (sampled_nodes.size() < oversample_fraction) {
							const count pos = sample(Aux::Random::getURNG());
							const node u = nodes_by_memberships[degree][pos];

							if (std::find(sampled_nodes.begin(), sampled_nodes.end(), u) == sampled_nodes.end()) {
								sampled_nodes.push_back(u);
							}
						}
					}

					assert(sampled_nodes.size() == nodes_with_fractional_slots.size() + 1);

					verifyInvariants();

					additionallyRemovedNode = sampled_nodes.front();
					removeNode(additionallyRemovedNode, false);

					verifyInvariants();

					// add a fraction to oversample_fraction nodes of nodes_of_degree
					for (count i = 0; i < nodes_with_fractional_slots.size(); ++i) {
						assert(nodes_with_fractional_slots[i] == none);

						const node u = sampled_nodes[i + 1];

						nodes_with_fractional_slots[i] = u;

						assert(nodes[u].degree == degree);

						for (count j = 0; j < u_additional_full_slots; ++j) {
							leaveCommunity(u);
						}

						addFraction(u, u_fractional_slot);
					}
				}
			}
			return additionallyRemovedNode;
		}

		node BucketSampling::removeNode(node nodeId) {
			assert(nodes[nodeId].degree > 0);
			if (nodes[nodeId].degree == 0) {
				throw std::runtime_error("Error, node has already been removed");
			}

			sumOfDesiredMemberships -= nodes[nodeId].degree;
			node result = removeNode(nodeId, true);
			verifyInvariants(true);
			return result;
		}

	}
}
