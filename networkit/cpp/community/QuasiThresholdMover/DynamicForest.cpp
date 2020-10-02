/*
 *
 */

#include <networkit/community/QuasiThresholdMover/DynamicForest.hpp>
#include <networkit/graph/GraphTools.hpp>

namespace NetworKit {
namespace QuasiThresholdMoving {

DynamicForest::DynamicForest() {}

// deprecated
DynamicForest::DynamicForest(const Graph &G, const std::vector<node> &parents)
    : path_membership(G.upperNodeIdBound(), none), path_pos(G.upperNodeIdBound(), 0), freeList(),
      paths(G.upperNodeIdBound(), SimplePath()) {
    Graph forest;
    if (parents.size() == G.upperNodeIdBound()) {
        forest = Graph(GraphTools::copyNodes(G), false, true);
        G.forNodes([&](node u) {
            if (parents[u] != none) {
                forest.addEdge(u, parents[u]);
            }
        });
    } else {
        forest = G;
    }
    // at first every node is in its own path
    std::iota(path_membership.begin(), path_membership.end(), 0);

    // build up parent/child relations
    forest.forNodes([&](node u) {
        paths[path(u)].pathNodes.push_back(u);
        if (forest.degreeOut(u) == 0) {
            paths[path(u)].posInParent = roots.size();
            roots.push_back(path(u));
            paths[path(u)].depth = 0;
        }
        forest.forEdgesOf(u, [&](node v) {
            paths[path(u)].parent = path(v);
            paths[path(u)].posInParent = paths[path(v)].childPaths.size();
            paths[path(v)].childPaths.push_back(path(u));
        });
    });
    // union nodes in simple paths by dfs
    pathDfsFrom(
        none,
        [&](pid subtreePath) {
            if (subtreePath != none) {
                pid p = paths[subtreePath].parent;
                if (p != none && paths[p].childPaths.size() == 1) {
                    unionPaths(p, subtreePath);
                }
            }
        },
        [](pid) {});
    updateDepthInSubtree(none);
    assert(pathsValid());
    TRACE("Dynamic Forest constructed");
}

DynamicForest::DynamicForest(const std::vector<node> &parents)
    : path_membership(parents.size(), none), path_pos(parents.size(), 0), freeList(),
      paths(parents.size(), SimplePath()) {
    // at first every node is in its own path
    std::iota(path_membership.begin(), path_membership.end(), 0);
    // build up parent/child relations
    for (node u = 0; u < parents.size(); u++) {
        paths[path(u)].pathNodes.push_back(u);
        node p = parents[u];
        if (p == none || p == u) {
            paths[path(u)].posInParent = roots.size();
            roots.push_back(path(u));
            paths[path(u)].depth = 0;
        } else {
            paths[path(u)].parent = path(p);
            paths[path(u)].posInParent = paths[path(p)].childPaths.size();
            paths[path(p)].childPaths.push_back(path(u));
        }
    }

    // union nodes in simple paths by dfs
    pathDfsFrom(
        none,
        [&](pid subtreePath) {
            if (subtreePath != none) {
                pid p = paths[subtreePath].parent;
                if (p != none && paths[p].childPaths.size() == 1) {
                    unionPaths(p, subtreePath);
                }
            }
        },
        [](pid) {});
    updateDepthInSubtree(none);
    assert(pathsValid());
    TRACE("Dynamic Forest constructed");
}

void DynamicForest::updateDepthInSubtree(pid start) {
    pathDfsFrom(
        start,
        [&](pid sp) {
            if (sp != none) {
                pid p = paths[sp].parent;
                if (p != none) {
                    paths[sp].depth = paths[p].depth + paths[p].length();
                } else {
                    paths[sp].depth = 0;
                }
            }
        },
        [](pid) {});
}

void DynamicForest::isolate(node u) {
    node oldParent = parent(u);
    if (paths[path(u)].length() == 1) {
        isolatePath(path(u));
    } else {
        isolateNode(u);
    }
    assert(pathsValid());
}

// isolate complete simple path
void DynamicForest::isolatePath(pid sp) {
    std::vector<pid> oldChildren = paths[sp].childPaths;
    assert(oldChildren.size() != 1);
    pid oldParent = paths[sp].parent;
    setParentPath(sp, none);
    // union paths if exactly one sibling is left back
    if (oldChildren.size() == 0 && oldParent != none && paths[oldParent].childPaths.size() == 1) {
        pid child = paths[oldParent].childPaths[0];
        unionPaths(oldParent, child);
        updateDepthInSubtree(child);
    } else {
        for (pid child : oldChildren) {
            setParentPath(child, oldParent);
            updateDepthInSubtree(child);
        }
    }
    updateDepthInSubtree(sp);

#ifndef NDEBUG
    for (pid i = 0; i < paths.size(); i++) {
        assert(paths[i].parent != sp);
        for (pid c : paths[i].childPaths) {
            assert(c != sp);
        }
    }
#endif
}

// isolate a node from a path of size >= 2
void DynamicForest::isolateNode(node u) {
    pid sp = path(u);
    index oldPos = path_pos[u];
    assert(paths[sp].length() >= 2);
    // update positions
    for (index i = oldPos + 1; i < paths[sp].length(); i++) {
        path_pos[paths[sp].pathNodes[i]] = i - 1;
    }
    paths[sp].pathNodes.erase(paths[sp].pathNodes.begin() + oldPos);
    // create new isolated path
    pid np = newPath();
    paths[np].posInParent = roots.size();
    roots.push_back(np);
    addToPath(u, np);
    updateDepthInSubtree(sp);
}

void DynamicForest::moveUpNeighbor(node neighbor, node referenceNode) {
    pid sp = path(neighbor);
    if (paths[sp].length() > 1) {
        TRACE("Move up ", neighbor);
        // update referenceNode if necessary
        if (paths[sp].referenceNode != referenceNode) {
            paths[sp].referenceNode = referenceNode;
            paths[sp].neighborCount = 0;
        }
        if (paths[sp].neighborCount >= paths[sp].length())
            return; // already all neighbors considered
        index oldPos = path_pos[neighbor];
        index neighborPos = paths[sp].length() - 1 - paths[sp].neighborCount;
        if (oldPos > neighborPos)
            return; // neighbor already considered
        paths[sp].neighborCount++;
        if (oldPos == neighborPos)
            return; // neighbor was not considered but is at right position
        node firstNonNeighbor = paths[sp].pathNodes[neighborPos];
        std::swap(paths[sp].pathNodes[path_pos[firstNonNeighbor]],
                  paths[sp].pathNodes[path_pos[neighbor]]);
        std::swap(path_pos[firstNonNeighbor], path_pos[neighbor]);
    }
    assert(pathsValid());
}

void DynamicForest::splitPath(pid sp, index splitPos) {
    if (paths[sp].length() <= 1 || splitPos == 0)
        return; // nothing to split
    assert(splitPos < paths[sp].length());
    pid oldParent = paths[sp].parent;
    index oldPosInParent = paths[sp].posInParent;
    count oldDepth = paths[sp].depth;

    pid np = newPath();
    // move upper nodes to new path
    for (int i = splitPos; i < paths[sp].length(); i++) {
        addToPath(paths[sp].pathNodes[i], np);
    }
    paths[sp].pathNodes.erase(paths[sp].pathNodes.begin() + splitPos, paths[sp].pathNodes.end());
    // update tree structure
    paths[sp].parent = np;
    paths[sp].posInParent = 0;
    paths[np].parent = oldParent;
    paths[np].posInParent = oldPosInParent;
    if (oldParent != none) {
        paths[oldParent].childPaths[oldPosInParent] = np;
    } else {
        roots[oldPosInParent] = np;
    }
    paths[np].childPaths.push_back(sp);
    paths[np].depth = oldDepth;
    paths[sp].depth = oldDepth + paths[np].length();
    assert(parent(paths[sp].upperEnd()) == paths[np].lowerEnd());
    assert(children(paths[np].lowerEnd()).size() == 1);
    assert(children(paths[np].lowerEnd())[0] == paths[sp].upperEnd());
}

void DynamicForest::unionPaths(pid upperPath, pid lowerPath) {
    if (upperPath == lowerPath) {
        return;
    }
    assert(paths[upperPath].childPaths.size() == 1);
    assert(paths[upperPath].childPaths[0] == lowerPath);

    pid upperParent = paths[upperPath].parent;
    pid upperPos = paths[upperPath].posInParent;

    // update tree structure
    paths[upperPath].childPaths.pop_back();
    paths[lowerPath].parent = upperParent;
    paths[lowerPath].posInParent = upperPos;
    if (upperParent != none) {
        paths[upperParent].childPaths[upperPos] = lowerPath;
    } else {
        roots[upperPos] = lowerPath;
    }
    // move nodes of upper path to lower path
    for (node u : paths[upperPath].pathNodes) {
        addToPath(u, lowerPath);
    }
    paths[upperPath].pathNodes.clear();
    deletePath(upperPath);
#ifndef NDEBUG
    if (paths[lowerPath].parent != none) {
        assert(paths[paths[lowerPath].parent].childPaths[upperPos] == lowerPath);
    } else {
        assert(roots[upperPos] == lowerPath);
    }
#endif
}

void DynamicForest::moveToPosition(node u, node p, const std::vector<node> &adoptedChildren) {
    // check that the node is isolated
    pid parentPath = path(p);
    pid oldPath = path(u);
#ifndef NDEBUG
    assert(paths[oldPath].parent == none);
    assert(parent(u) == none);
    assert(childCount(u) == 0);
    assert(paths[path(u)].length() == 1);
    // check that all children are adopted ones
    if (p != none) {
#ifndef NDEBUG
        std::vector<node> oldChildren = children(p);
        for (node c : adoptedChildren) {
            assert(std::find(oldChildren.begin(), oldChildren.end(), c) != oldChildren.end());
        }
#endif
    } else {
        for (node c : adoptedChildren) {
            assert(std::find(roots.begin(), roots.end(), path(c)) != roots.end());
        }
    }
#endif
    // if all children are adopted, position in path does not matter, insert on top
    if (adoptedChildren.size() == childCount(p)) {
        path_pos[u] = paths[parentPath].pathNodes.size();
        paths[parentPath].pathNodes.push_back(u);
        path_membership[u] = parentPath;

        // update tree structure
        index oldTreePos = paths[oldPath].posInParent;
        roots[oldTreePos] = roots.back();
        roots.pop_back();
        paths[roots[oldTreePos]].posInParent = oldTreePos;

        paths[oldPath].pathNodes.pop_back();
        deletePath(oldPath);
    } else {
        // if at least one child is not adopted, parent needs to be lower end
        if (p != none && !isLowerEnd(p)) {
            splitPath(parentPath, path_pos[p]);
            parentPath = path(p);
        }
        // if exactly one child is adopted, node can just be added to that path
        if (adoptedChildren.size() == 1) {
            // update Tree structure
            index oldTreePos = paths[oldPath].posInParent;
            roots[oldTreePos] = roots.back();
            roots.pop_back();
            paths[roots[oldTreePos]].posInParent = oldTreePos;

            paths[oldPath].pathNodes.pop_back();
            deletePath(oldPath);

            addToPath(u, path(adoptedChildren[0]));
            // otherwise, insert in simple-path tree structure
        } else {
            setParentPath(oldPath, parentPath);
            for (node child : adoptedChildren) {
                setParentPath(path(child), path(u));
            }
        }
    }
    updateDepthInSubtree(path(u));
    assert(pathsValid());
}

bool DynamicForest::pathsValid() {
    // check that parent/child relations for paths are valid
    for (pid sp = 0; sp < paths.size(); sp++) {
        if (paths[sp].pathNodes.size() > 0) {
            const std::vector<pid> &cps = paths[sp].childPaths;
            for (index i = 0; i < cps.size(); i++) {
                pid child = cps[i];
                assert(paths[child].parent == sp);
                assert(paths[child].posInParent == i);
            }
        }
    }
    for (index i = 0; i < roots.size(); i++) {
        pid rootPath = roots[i];
        assert(paths[rootPath].parent == none);
        assert(paths[rootPath].posInParent == i);
    }

    for (node u = 0; u < path_membership.size(); u++) {
        const std::vector<node> childNodes = children(u);
        assert(childNodes.size() == childCount(u));
        // check that parent/child realtionships for nodes are proper
        for (index i = 0; i < childCount(u); i++) {
            node c = childNodes[i];
            assert(u != c);
            if (i > 0) {
                assert(path(u) != path(c));
            }
            assert(posInParent(c) == i);
            assert(parent(c) == u);
        }
        // check that  path_membership and positions for nodes are proper
        pid sp = path(u);
        if (paths[sp].pathNodes.size() > 0) {
            assert(paths[sp].pathNodes[path_pos[u]] == u);
            for (index i = 0; i < paths[sp].pathNodes.size(); i++) {
                node v = paths[sp].pathNodes[i];
                assert(path_pos[v] == i);
                assert(path(v) == sp);
            }
            const std::vector<pid> &cps = paths[sp].childPaths;
            for (pid i = 0; i < paths.size(); i++) {
                if (paths[i].parent == sp) {
                    assert(std::find(cps.begin(), cps.end(), i) != cps.end());
                }
            }
        }
    }
    // check that depth is proper for paths
    pathDfsFrom(
        none,
        [&](pid sp) {
            if (sp != none) {
                pid p = paths[sp].parent;
                if (p != none) {
                    if (paths[sp].depth != paths[p].depth + paths[p].length()) {
                    }
                    assert(paths[sp].depth == paths[p].depth + paths[p].length());
                } else {
                    if (paths[sp].depth != 0) {
                    }
                    assert(paths[sp].depth == 0);
                }
            }
        },
        [](pid) {});

    // check that depth is proper for nodes
    dfsFrom(
        none,
        [&](node c) {
            if (c != none && parent(c) != none) {
                assert(depth(c) == depth(parent(c)) + 1);
                if (childCount(c) == 1) {
                    assert(path_membership[c] == path_membership[children(c)[0]]);
                }
            }
        },
        [](node) {});

    // check that free list is proper
    for (pid freePlace : freeList) {
        assert(paths[freePlace].parent == none);
        assert(paths[freePlace].referenceNode == none);
        assert(paths[freePlace].pathNodes.size() == 0);
        assert(paths[freePlace].childPaths.size() == 0);
        assert(paths[freePlace].neighborCount == 0);
        assert(paths[freePlace].depth == 0);
        assert(paths[freePlace].posInParent == 0);
    }
    return 1;
}

std::string DynamicForest::printPaths() const {
    std::stringstream ss;
    for (node u = 0; u < path_membership.size(); u++) {
        ss << "{" << u << " ";
        ss << "[";
        std::vector<node> pNodes = paths[path(u)].pathNodes;
        for (node u : pNodes) {
            ss << u << " ";
        }
        ss << "]";
        if (path_pos[u] != none) {
            ss << " pos " << path_pos[u];
        }
        ss << "} ";
    }
    return ss.str();
}

Graph DynamicForest::toGraph() const {
    Graph result(path_membership.size(), false, true);

    for (pid r : roots) {
        dfsFrom(
            paths[r].upperEnd(), [](node) {},
            [&](node u) {
                if (parent(u) != none) {
                    assert(u < path_membership.size());
                    assert(parent(u) < path_membership.size());
                    result.addEdge(u, parent(u));
                }
            });
    }

    return result;
}

void DynamicForest::deletePath(DynamicForest::pid i) {
// check that path got isolated from tree structure
#ifndef NDEBUG
    assert(i != none);
    assert(paths[i].childPaths.size() == 0);
    for (pid sp = 0; sp < paths.size(); sp++) {
        assert(paths[sp].parent != i);
        for (pid c : paths[sp].childPaths) {
            assert(c != i);
        }
    }
    assert(std::find(roots.begin(), roots.end(), i) == roots.end());
#endif
    paths[i].reset();
    freeList.push_back(i);
}

DynamicForest::pid DynamicForest::newPath() {
    assert(!freeList.empty());
    pid freePlace = freeList.back();
    freeList.pop_back();
    return freePlace;
}

count DynamicForest::childCount(node u) const {
    if (u == none) {
        return roots.size();
    }
    if (!isLowerEnd(u)) {
        return 1;
    } else {
        return paths[path(u)].childPaths.size();
    }
}

index DynamicForest::posInParent(node u) const {
    if (u == none) {
        return none;
    } else if (!isUpperEnd(u)) {
        return 0;
    } else {
        return paths[path(u)].posInParent;
    }
}

} // namespace QuasiThresholdMoving
} // namespace NetworKit
