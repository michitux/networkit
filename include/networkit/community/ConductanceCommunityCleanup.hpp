// networkit-format
#ifndef NETWORKIT_COMMUNITY_CONDUCTANCE_COMMUNITY_CLEANUP_HPP_
#define NETWORKIT_COMMUNITY_CONDUCTANCE_COMMUNITY_CLEANUP_HPP_

#include <tlx/container/d_ary_addressable_int_heap.hpp>

#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/structures/GlobalCommunity.hpp>

namespace NetworKit {
class ConductanceCommunityCleanup final : public Algorithm {
public:
    ConductanceCommunityCleanup(const Graph &g);

    ConductanceCommunityCleanup(const ConductanceCommunityCleanup &other);

    template <typename InputIt>
    void setCommunity(InputIt begin, InputIt end) {
        community.clear();
        community.init(begin, end);
        hasRun = false;
    }

    void run() override;

    template <typename F>
    void forCommunityNodes(F &&f) const {
        assureFinished();
        community.forCommunityNodes(f);
    }

    ConductanceCommunityCleanup &operator=(const ConductanceCommunityCleanup &other);

    size_t size() const { return community.getSize(); }

private:
    double communityStrength(node u) const {
        return community.getInternalStrength(u) / weightedDegree[u];
    }

    class StrengthLess {
    public:
        StrengthLess(const ConductanceCommunityCleanup &cc) : cc(&cc) {}
        bool operator()(node a, node b) const {
            return cc->communityStrength(a) < cc->communityStrength(b);
        }

    private:
        const ConductanceCommunityCleanup *cc;
    };

    class StrengthGreater {
    public:
        StrengthGreater(const ConductanceCommunityCleanup &cc) : cc(&cc) {}
        bool operator()(node a, node b) const {
            return cc->communityStrength(a) > cc->communityStrength(b);
        }

    private:
        const ConductanceCommunityCleanup *cc;
    };

    const Graph *g;
    GlobalCommunity community;
    double totalVolume;
    std::vector<double> weightedDegree;

    tlx::DAryAddressableIntHeap<node, 4, StrengthLess> communityHeap;
    tlx::DAryAddressableIntHeap<node, 4, StrengthGreater> shellHeap;
};
} // namespace NetworKit

#endif // NETWORKIT_COMMUNITY_CONDUCTANCE_COMMUNITY_CLEANUP_HPP_
