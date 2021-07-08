// networkit-format

#include <networkit/community/ConductanceCommunityCleanup.hpp>

namespace NetworKit {
ConductanceCommunityCleanup::ConductanceCommunityCleanup(const Graph &g)
    : g(&g), community(g), totalVolume(0), weightedDegree(g.upperNodeIdBound()),
      communityHeap(StrengthLess(*this)), shellHeap(StrengthGreater(*this)) {
    totalVolume = g.totalEdgeWeight() * 2;

    g.parallelForNodes([&](node u) { weightedDegree[u] = g.weightedDegree(u, true); });
}

// This should be a standard copy constructor, but we do not copy
// the heaps as they contain a pointer to this.
ConductanceCommunityCleanup::ConductanceCommunityCleanup(const ConductanceCommunityCleanup &other)
    : Algorithm(other), g(other.g), community(other.community), totalVolume(other.totalVolume),
      weightedDegree(other.weightedDegree), communityHeap(StrengthLess(*this)),
      shellHeap(StrengthGreater(*this)) {}

ConductanceCommunityCleanup &
ConductanceCommunityCleanup::operator=(const ConductanceCommunityCleanup &other) {
    Algorithm::operator=(other);
    g = other.g;
    community = other.community;
    totalVolume = other.totalVolume;
    weightedDegree = other.weightedDegree;

    // shellHeap and communityHeap do not depend on the graph and should
    // always be in a cleared state. However, they do contain a pointer
    // to this, so we must not copy them.
    communityHeap.clear();
    shellHeap.clear();

    return *this;
}

void ConductanceCommunityCleanup::run() {
    size_t initialSize = community.getSize();

    TRACE("Starting cleanup of community of size ", initialSize);

    community.forNodesWithInternalStrength([&](node u, double) {
        if (community.contains(u)) {
            communityHeap.push(u);
        } else {
            shellHeap.push(u);
        }
    });

    auto conductance = [&](double cut, double volume) {
        return cut / std::min(volume, totalVolume - volume);
    };

    double currentConductance = conductance(community.getCut(), community.getVolume());

    count numAdded = 0, numRemoved = 0;

    auto tooManyChanged = [&]() -> bool { return (numAdded + numRemoved) > initialSize; };

    while (!(communityHeap.empty() && shellHeap.empty()) && !tooManyChanged()) {
        double cutAdd = 0, volAdd = 0, cutRemove = 0, volRemove = 0;
        double conductanceAdd = std::numeric_limits<double>::max();
        double conductanceRemove = std::numeric_limits<double>::max();

        if (!communityHeap.empty()) {
            node u = communityHeap.top();
            // FIXME: this is not correct for graphs with self-loops!
            cutRemove =
                community.getCut() + 2 * community.getInternalStrength(u) - weightedDegree[u];
            volRemove = community.getVolume() - weightedDegree[u];

            conductanceRemove = conductance(cutRemove, volRemove);
        }

        if (!shellHeap.empty()) {
            node u = shellHeap.top();
            cutAdd = community.getCut() - 2 * community.getInternalStrength(u) + weightedDegree[u];
            volAdd = community.getVolume() + weightedDegree[u];

            conductanceAdd = conductance(cutAdd, volAdd);
        }

        if (conductanceRemove <= conductanceAdd) {
            node u = communityHeap.top();
            communityHeap.pop();

            if (conductanceRemove < currentConductance) {
                currentConductance = conductanceRemove;

                community.removeNode(u, [&](node v, double internalStrength) {
                    if (community.contains(v)) {
                        // This may insert the node again if it is not in the queue.
                        // This is intentional here as the node just lost a neighbor
                        communityHeap.update(v);
                    } else if (shellHeap.contains(v)) {
                        if (internalStrength > 0) {
                            shellHeap.update(v);
                        } else {
                            shellHeap.remove(v);
                        }
                    }
                });

                ++numRemoved;

                assert(currentConductance
                       == conductance(community.getCut(), community.getVolume()));
            }
        } else {
            node u = shellHeap.top();
            shellHeap.pop();

            assert(!community.contains(u));

            if (conductanceAdd < currentConductance) {
                currentConductance = conductanceAdd;
                ++numAdded;

                community.addNode(u, [&](node v, double) {
                    // Update heaps
                    if (!community.contains(v)) {
                        // This may insert the node if it is not in the queue
                        // which is intential here as this node got a higher
                        // score now
                        shellHeap.update(v);
                    } else if (communityHeap.contains(v)) {
                        // Only update the community heap if the node is already
                        // in it, as this node has a new internal neighbor it
                        // seems unlikely that it should be removed now
                        communityHeap.update(v);
                    }
                });

                assert(currentConductance
                       == conductance(community.getCut(), community.getVolume()));
            }
        }
    }

    // FIXME: clear is linear in the key space!
    while (!communityHeap.empty()) {
        communityHeap.pop();
    }
    while (!shellHeap.empty()) {
        shellHeap.pop();
    }

    if (tooManyChanged()) {
        community.clear();
    }

    hasRun = true;
}

} // namespace NetworKit
