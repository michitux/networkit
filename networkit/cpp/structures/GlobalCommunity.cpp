// networkit-format

#include <limits>

#include <networkit/structures/GlobalCommunity.hpp>

namespace NetworKit {

GlobalCommunity::GlobalCommunity(const Graph &g)
    : g(&g), inCommunity(g.upperNodeIdBound(), CommunityStatus::unused),
      internalStrength(g.upperNodeIdBound(), std::numeric_limits<double>::infinity()), volume(0),
      cut(0), size(0) {}

void GlobalCommunity::clear() {
    inCommunity.reset();
    internalStrength.reset();
    volume = 0;
    cut = 0;
    size = 0;
}

void GlobalCommunity::checkConsistency() const {
#ifndef NDEBUG
#ifdef NETWORKIT_SANITY_CHECKS
    double volumeDebug = 0, cutDebug = 0;
    size_t sizeDebug = 0;

    g->forNodes([&](node u) {
        double debugInternalStrength = 0;

        g->forNeighborsOf(u, [&](node v, edgeweight ew) {
            if (u == v)
                return;

            if (inCommunity[v] == CommunityStatus::inside) {
                debugInternalStrength += ew;
            }
        });

        if (internalStrength.indexIsUsed(u)) {
            assert(debugInternalStrength == internalStrength[u]);
        } else {
            assert(debugInternalStrength == 0);
        }

        if (inCommunity[u] == CommunityStatus::inside) {
            ++sizeDebug;
            g->forNeighborsOf(u, [&](node v, edgeweight ew) {
                volumeDebug += ew;

                if (u == v) {
                    volumeDebug += ew;
                    return;
                }

                if (inCommunity[v] != CommunityStatus::inside) {
                    cutDebug += ew;
                }
            });
        }
    });

    assert(volumeDebug == volume);
    assert(cutDebug == cut);
    assert(sizeDebug == size);
#endif
#endif
}
} // namespace NetworKit
