// networkit-format
#ifndef NETWORKIT_STRUCTURES_GLOBAL_COMMUNITY_HPP_
#define NETWORKIT_STRUCTURES_GLOBAL_COMMUNITY_HPP_

#include <vector>

#include <networkit/auxiliary/SparseVector.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {

/**
 * Maintain a single community that can be expanded and shrinked
 *
 * In contrast to LocalCommunity, this data structure uses memory
 * linear in the number of nodes of the graph to be more efficient.
 * This is in particular beneficial if several communities re-use the
 * same data structure, as then the initialization overhead that is
 * linear in the number of nodes occurs only once.
 *
 * The feature set is currently limited due to limited use cases but
 * may be expanded in the future.
 */
class GlobalCommunity final {
public:
    /**
     * Initialize the global community data structure for the given graph.
     */
    GlobalCommunity(const Graph &g);

    /**
     * Initialize the data structure for a community
     *
     * @param begin Iterator pointing to the first node to insert
     * @param end Past the end iterator of the nodes to insert
     */
    template <typename InputIt>
    void init(InputIt begin, InputIt end) {
        assert(inCommunity.size() == 0);
        assert(internalStrength.size() == 0);

        for (InputIt it = begin; it != end; ++it) {
            inCommunity.insertOrSet(*it, CommunityStatus::inside);
            ++size;
        }

        forCommunityNodes([&](node u) {
            g->forNeighborsOf(u, [&](node v, edgeweight ew) {
                volume += ew;

                if (u == v) {
                    volume += ew;
                    return;
                }

                if (!internalStrength.indexIsUsed(v)) {
                    internalStrength.insert(v, ew);
                } else {
                    internalStrength[v] += ew;
                }

                if (inCommunity[v] != CommunityStatus::inside) {
                    cut += ew;
                }
            });
        });

        checkConsistency();
    }

    /**
     * Add a single node to the community
     *
     * @param u The node to add
     * @param onNeighborChange Callback to call for every neighbor of
     * the node after its internal strength has been updated
     */
    template <typename F>
    void addNode(node u, F &&onNeighborChange) {
        assert(inCommunity[u] != CommunityStatus::inside);

        inCommunity.insertOrSet(u, CommunityStatus::inside);

        ++size;

        g->forNeighborsOf(u, [&](node v, edgeweight ew) {
            volume += ew;

            if (u == v) {
                volume += ew;
                return;
            }

            if (!internalStrength.indexIsUsed(v)) {
                internalStrength.insert(v, ew);
            } else {
                internalStrength[v] += ew;
            }

            if (inCommunity[v] == CommunityStatus::inside) {
                cut -= ew;
            } else {
                cut += ew;
            }

            onNeighborChange(v, internalStrength[v]);
        });
    }

    /**
     * Remove a single node from the community
     *
     * @param u The node to remove
     * @param onNeighborChange Callback to call for every neighbor of
     * the node after its internal strength has been updated
     */
    template <typename F>
    void removeNode(node u, F &&onNeighborChange) {
        assert(inCommunity[u] == CommunityStatus::inside);

        inCommunity[u] = CommunityStatus::removed;

        --size;

        g->forNeighborsOf(u, [&](node v, edgeweight ew) {
            volume -= ew;

            if (u == v) {
                volume -= ew;
                return;
            }

            internalStrength[v] -= ew;

            if (inCommunity[v] == CommunityStatus::inside) {
                cut += ew;
            } else {
                cut -= ew;
            }

            onNeighborChange(v, internalStrength[v]);
        });
    }

    /**
     * Iterate over all nodes that have an internal strength set
     *
     * Note that if nodes have been removed, the callback may be
     * called for nodes with strength 0.
     *
     * @param f The callback to call for every node with the node and the internal strength
     */
    template <typename F>
    void forNodesWithInternalStrength(F &&f) const {
        internalStrength.forElements(std::forward<F>(f));
    }

    /**
     * Iterate over all nodes in the community
     *
     * Note that currently the running time is linear in the total number of number of nodes that
     * were part of the community at some point after the last time clear has been called.
     *
     * @param f The callback to call for every node in the community
     */
    template <typename F>
    void forCommunityNodes(F &&f) const {
        inCommunity.forElements([&](node u, CommunityStatus status) {
            if (status == CommunityStatus::inside) {
                f(u);
            }
        });
    }

    /**
     * Query if a node is contained in the community
     *
     * @param u The node to check
     * @return If @a u is in the community
     */

    bool contains(node u) const { return inCommunity[u] == CommunityStatus::inside; }

    /**
     * Get the internal strength of a node
     *
     * The value is undefined if the node is neither in the community or a neighbor of the
     * community.
     *
     * @param u The node to check
     * @return The internal strength of @a u
     */
    double getInternalStrength(node u) const { return internalStrength[u]; }

    /**
     * Reset all internal data structures such that a new community can be managed
     */
    void clear();

    /**
     * Get the cut of the community
     *
     * @return The cut of the community
     */
    double getCut() const { return cut; }

    /**
     * Get the volume of the community
     *
     * @return The volume of the community
     */
    double getVolume() const { return volume; }

    /**
     * Get the number of nodes in the community
     *
     * @return Thes size of the community
     */
    size_t getSize() const { return size; }

private:
    /**
     * Check if all data structures and calculated values are consistent
     */
    void checkConsistency() const;

    enum class CommunityStatus : uint8_t { unused = 0, inside = 1, removed = 2 };

    const Graph *g;
    SparseVector<CommunityStatus> inCommunity;
    SparseVector<double> internalStrength;

    double volume;
    double cut;
    size_t size;
};

} // namespace NetworKit

#endif // NETWORKIT_STRUCTURES_GLOBAL_COMMUNITY_HPP_
