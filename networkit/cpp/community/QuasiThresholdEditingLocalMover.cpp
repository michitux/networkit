/*
 *
 */

#include "QuasiThresholdEditingLocalMover.h"
#include "../generators/TreeReachabilityGraphGenerator.h"
#include <unordered_set>
#include "../graph/DynamicForest.h"

NetworKit::QuasiThresholdEditingLocalMover::QuasiThresholdEditingLocalMover(const NetworKit::Graph &G, const std::vector< NetworKit::node > &parent, NetworKit::count maxIterations)
: G(G), maxIterations(maxIterations), numEdits(none) {
	forest = Graph(G.copyNodes(), false, true);
	G.forNodes([&](node u) {
		if (parent[u] != none) {
			forest.addEdge(u, parent[u]);
		}
	});
}

void NetworKit::QuasiThresholdEditingLocalMover::run() {
	/*
	Idea: for each node check if it can be inserted below any other node such that the number of edits is reduced.
	Additionally check if the node should become a root node and add any number of previous roots as children.
	For each child of the new parent check if it can be inserted below the moved node. Question: is this possible such that the total time for the move is linear?

	For a node that shall be moved, mark all neighbors in the original graph.
	For just checking if a node can be inserted as leaf under an anchor node or as parent of an anchor node only the neighbors of that anchor node need to be iterated
	and checked for the mark, then the number of necessary edits can be determined.

	Problem: In order to keep these checks easy we should remove the node from the graph. Therefore insertion in the previous position must remain possible. This means that
	we need to be able to add arbitrary children the node previously had. We find them when iterating over the neighbors of the parent.
	But can we determine which of them should be added?

	Another idea: Keep an efficient forest data structure (or simply rooted trees in a graph + a list of roots?) that allows walking up and down in the tree.
	Then do DFS in order to discover possible inserts/deletes.
	This should actually allow to determine which children to take. Furthermore then probably we don't need to store the inserted edges explicitly but we can determine
	them again whenever we delete a node as we iterate over its neighbors anyway and we can in the same turn also iterate over the tree and determine which neighbors were there
	and which were not by checking for marked nodes.
	*/

	DynamicForest dynamicForest(forest);

	// count number of edits that are needed with the initial given forest
	count numExistingEdges = 0;
	count numMissingEdges = 0;
	std::vector<bool> marker(G.upperNodeIdBound());

	dynamicForest.forChildrenOf(none, [&](node r) {
		count depth = 0;
		dynamicForest.dfsFrom(r,
		[&](node u) { // on enter
			count upperNeighbors = 0;

			G.forNeighborsOf(u, [&](node v) {
				if (marker[v]) {
					++upperNeighbors;
				}
			});

			numExistingEdges += upperNeighbors;
			numMissingEdges += depth - upperNeighbors;

			marker[u] = true;
			depth += 1;
		},
		[&](node u) { // on exit
			marker[u] = false;
			depth -= 1;
		});
	});

	numEdits = numMissingEdges + (G.numberOfEdges() - numExistingEdges);

	bool hasMoved = true;

	std::vector<count> numNeighbors(G.upperNodeIdBound(), 0);

	for (count i = 0; hasMoved && i < maxIterations; ++i) {
		hasMoved = false;

		G.forNodesInRandomOrder([&](node nodeToMove) {
			G.forEdgesOf(nodeToMove, [&](node v) {
				marker[v] = true;
			});

			std::vector<count> missingBelow(G.upperNodeIdBound(), 0), missingAbove(G.upperNodeIdBound(), 0), existingBelow(G.upperNodeIdBound(), 0), existingAbove(G.upperNodeIdBound(), 0);

			dynamicForest.forChildrenOf(none, [&](node r) {
				dynamicForest.dfsFrom(r,
				[&](node u) {
					if (u != nodeToMove) {
						missingBelow[u] = missingAbove[u] = 1 - marker[u];
						existingBelow[u] = existingAbove[u] = marker[u];
					}
					node p = dynamicForest.parent(u);
					if (p != none) {
						missingAbove[u] += missingAbove[p];
						existingAbove[u] += existingAbove[p];
					}
				},
				[&](node u) {
					node p = dynamicForest.parent(u);
					if (p != none) {
						missingBelow[p] += missingBelow[u];
						existingBelow[p] += existingBelow[u];
					}
				});
			});

			count curEdits = G.degree(nodeToMove) - existingAbove[nodeToMove] - existingBelow[nodeToMove] + missingAbove[nodeToMove] + missingBelow[nodeToMove];
			count minEdits = curEdits;
			std::vector<node> curChildren(dynamicForest.children(nodeToMove));
			std::vector<node> minChildren;
			node curParent = dynamicForest.parent(nodeToMove);
			node minParent = curParent;

			dynamicForest.isolate(nodeToMove);

			auto tryEditBelow = [&](node p) {
				if (p == nodeToMove) return;

				count edits = G.degree(nodeToMove);
				if (p != none) {
					edits += missingAbove[p];
					edits -= existingAbove[p];
				}

				std::vector<node> children;
				dynamicForest.forChildrenOf(p, [&](node c) {
					if (c == nodeToMove) return;
					if (existingBelow[c] > missingBelow[c]) { // TODO try >= (more children...)
						children.emplace_back(c);
						edits -= existingBelow[c] - missingBelow[c];
					}
				});

				if (edits < minEdits) {
					minEdits = edits;
					minChildren = std::move(children);
					minParent = p;
				}

			};

			G.forNodes(tryEditBelow);
			tryEditBelow(none);

			dynamicForest.setParent(nodeToMove, curParent);
			for (node c : curChildren) {
				dynamicForest.setParent(c, nodeToMove);
			}

			G.forEdgesOf(nodeToMove, [&](node v) {
				marker[v] = false;
			});


			count subtreeSize = 0;
			dynamicForest.dfsFrom(nodeToMove, [&](node d) {
				marker[d] = true;
				++subtreeSize;
			}, [](node d) {});

			count subtreeExtDegree = 0;
			dynamicForest.dfsFrom(nodeToMove, [&](node d) {
				G.forNeighborsOf(d, [&](node v) {
					if (!marker[v]) {
						++numNeighbors[v];
						++subtreeExtDegree;
					}
				});
			}, [](node d) {});

			dynamicForest.forChildrenOf(none, [&](node r) {
				dynamicForest.dfsFrom(r,
				[&](node u) {
					if (!marker[u]) {
						missingAbove[u] = subtreeSize - numNeighbors[u];
						existingAbove[u] = numNeighbors[u];
						node p = dynamicForest.parent(u);
						if (p != none) {
							missingAbove[u] += missingAbove[p];
							existingAbove[u] += existingAbove[p];
						}
					}
				},
				[](node u) {
				});
			});

			// virtually remove the subtree of u from the forest concerning the missing/existingBelow-counters
			for (node u = curParent; u != none; u = dynamicForest.parent(u)) {
				existingBelow[u] -= existingBelow[nodeToMove];
				missingBelow[u] -= missingBelow[nodeToMove];
			}

			bool moveWithSubtree = false;

			// calculate how many edits the whole subtree currently needs
			count curSubtreeEdits = subtreeExtDegree;
			if (curParent != none) { // FIXME here we ignore edits inside the tree as they do not change. Is this okay?
				curSubtreeEdits += missingAbove[curParent];
				curSubtreeEdits -= existingAbove[curParent];
			}

			// calculate the number of saved edits as comparing the absolute number of edits doesn't make sense
			count savedEdits = curEdits - minEdits;

			auto trySubtreeEditBelow = [&](node p) {
				if (p != none && marker[p]) return;

				count edits = subtreeExtDegree;
				if (p != none) {
					edits += missingAbove[p];
					edits -= existingAbove[p];
				}

				std::vector<node> children;
				dynamicForest.forChildrenOf(p, [&](node c) {
					if (marker[c]) return;
					if (existingBelow[c] >= missingBelow[c]) { // TODO try >= (more children...)
						children.emplace_back(c);
						edits -= existingBelow[c] - missingBelow[c];
					}
				});

				if (edits < curSubtreeEdits && savedEdits < curSubtreeEdits - edits) {
					minEdits = edits;
					minChildren = std::move(children);
					minParent = p;
					moveWithSubtree = true;
					savedEdits = curSubtreeEdits - edits;
				}

			};

			G.forNodes(trySubtreeEditBelow);
			trySubtreeEditBelow(none);

			dynamicForest.dfsFrom(nodeToMove, [&](node d) {
				marker[d] = false;
				G.forNeighborsOf(d, [&](node v) {
					numNeighbors[v] = 0;
				});
			}, [](node d) {});


			if (savedEdits > 0) {
				if (!moveWithSubtree) {
					dynamicForest.isolate(nodeToMove);
				}
				dynamicForest.setParent(nodeToMove, minParent);
				for (node c : minChildren) {
					dynamicForest.setParent(c, nodeToMove);
				}

				hasMoved = true;
				numEdits -= savedEdits;
			}

		});
	}

	forest = dynamicForest.toGraph();
}

NetworKit::count NetworKit::QuasiThresholdEditingLocalMover::getNumberOfEdits() {
	return numEdits;
}

NetworKit::Graph NetworKit::QuasiThresholdEditingLocalMover::getQuasiThresholdGraph() {
	TreeReachabilityGraphGenerator gen(forest);
	gen.run();
	return gen.getGraph();
}



