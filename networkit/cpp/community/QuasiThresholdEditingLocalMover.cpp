/*
 *
 */

#include "QuasiThresholdEditingLocalMover.h"
#include "../generators/TreeReachabilityGraphGenerator.h"
#include <unordered_set>
#include "../graph/DynamicForest.h"

NetworKit::QuasiThresholdEditingLocalMover::QuasiThresholdEditingLocalMover(const NetworKit::Graph &G, const std::vector< NetworKit::node > &parent, NetworKit::count maxIterations, bool moveSubtrees)
: G(G), maxIterations(maxIterations), moveSubtrees(moveSubtrees), numEdits(none) {
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

	numEdits = countNumberOfEdits();
	usedIterations = 0;

	bool hasMoved = true;
	std::vector<bool> marker(G.upperNodeIdBound());

	std::vector<count> numNeighbors;

	if (moveSubtrees) {
		numNeighbors.resize(G.upperNodeIdBound(), 0);
	}

	std::vector<count> depth(G.upperNodeIdBound(), 0);
	dynamicForest.dfsFrom(none, [&](node c) {
		if (c != none && dynamicForest.parent(c) != none) {
			depth[c] = depth[dynamicForest.parent(c)] + 1;
		}
	}, [](node){});

	std::vector<node> neighborQueue, currentLevel, nextLevel, touchedNodes, lastVisitedDFSNode(G.upperNodeIdBound(), none), bestParentBelow(G.upperNodeIdBound(), none);
	G.parallelForNodes([&](node u) {
		lastVisitedDFSNode[u] = u;
	});
	std::vector<count> maxGain(G.upperNodeIdBound(), 0), editDifference(G.upperNodeIdBound(), 0);
	std::vector<bool> nodeTouched(G.upperNodeIdBound(), false);

	for (count i = 0; hasMoved && i < maxIterations; ++i) {
		hasMoved = false;

		G.forNodesInRandomOrder([&](node nodeToMove) {
			G.forEdgesOf(nodeToMove, [&](node v) {
				marker[v] = true;
			});

			// remove the node from its tree but store the old position.
			std::vector<node> curChildren(dynamicForest.children(nodeToMove));
			node curParent = dynamicForest.parent(nodeToMove);
			count curEdits = G.degree(nodeToMove);

			dynamicForest.dfsFrom(nodeToMove, [&](node c) {
				--depth[c];
				if (c != nodeToMove) {
					curEdits += 1 - 2 * marker[c];
				}
			}, [](node){});

			for (node p = dynamicForest.parent(nodeToMove); p != none; p = dynamicForest.parent(p)) {
				curEdits += 1 - 2 * marker[p];
			}

			dynamicForest.isolate(nodeToMove);

			depth[nodeToMove] = 0;

			G.forEdgesOf(nodeToMove, [&](node v) {
				neighborQueue.emplace_back(v);
			});

			std::stable_sort(neighborQueue.begin(), neighborQueue.end(), [&](node u, node v) {return depth[u] < depth[v];}); // the queue shall be used from the end

			count level = 0;
			if (!neighborQueue.empty()) {
				level = depth[neighborQueue.back()];
			}

			count bestParent = none;
			count rootMaxGain = 0, rootEdits = 0;

			auto processNode = [&](node u) {
				TRACE("Processing node ", u, " of depth ", depth[u], " (node to move: ", nodeToMove, ")");
				TRACE("Parent: ", dynamicForest.parent(u), ", children: ", dynamicForest.children(u));
				assert(u != nodeToMove);

				if (!nodeTouched[u]) {
					nodeTouched[u] = true;
					touchedNodes.emplace_back(u);
					assert(marker[u]); // only marked neighbors may be processed without having been touched before
				}

				count sumPositiveEdits = editDifference[u];

				assert(editDifference[u] != none);

				editDifference[u] += marker[u];
				editDifference[u] -= 1 - marker[u]; // if (marker[u]) { ++editDifference[u]; } else { --editDifference[u]; }

				TRACE("Edit difference before descending: ", editDifference[u]);

				assert(!marker[u] || editDifference[u] > 0);

				if (editDifference[u] != none) {
					assert(lastVisitedDFSNode[u] == u);

					node c = dynamicForest.nextDFSNodeOnEnter(u, u);

					while (c != u) {
						assert(depth[c] > level);

						if (!nodeTouched[c] || editDifference[c] == none) {
							--editDifference[u];

							// advance to the next starting point for the DFS search.
							c = lastVisitedDFSNode[c];

							if (editDifference[u] == none) {
								lastVisitedDFSNode[u] = c;
								break;
							}

							c = dynamicForest.nextDFSNodeOnEnter(c, u);
						} else {
							node p = dynamicForest.parent(c);
							c = dynamicForest.nextChild(c, p);

							while (c == p && c != u) {
								p = dynamicForest.parent(p);
								c = dynamicForest.nextChild(c, p);
							}
						}
					}
				}

				TRACE("Edit difference after descending: ", editDifference[u]);

				if (sumPositiveEdits > maxGain[u] || maxGain[u] == 0) {
					maxGain[u] = sumPositiveEdits;
					bestParentBelow[u] = u;
				}

				assert(maxGain[u] != none);

				maxGain[u] += marker[u];

				if (maxGain[u] > 0) {
					maxGain[u] -= 1 - marker[u];
				}

				TRACE("Maximum gain at ", u, ": ", maxGain[u]);

				if (maxGain[u] > 0 || (editDifference[u] != none && editDifference[u] != 0)) {
					node p = dynamicForest.parent(u);

					if (p != none) {
						if (!nodeTouched[p]) {
							nodeTouched[p] = true;
							touchedNodes.push_back(p);
							nextLevel.push_back(p);
						}

						if (editDifference[u] != none) {
							assert(editDifference[u] <= maxGain[u]);
							editDifference[p] += editDifference[u];
						}

						if (maxGain[u] > maxGain[p]) {
							maxGain[p] = maxGain[u];
							bestParentBelow[p] = bestParentBelow[u];
						}
					} else {
						if (maxGain[u] > rootMaxGain) {
							rootMaxGain = maxGain[u];
							bestParent = bestParentBelow[u];
						}

						if (editDifference[u] != none) {
							rootEdits += editDifference[u];
						}
					}
				}

				if (dynamicForest.children(u).empty()) { assert(editDifference[u] == 1); }
			};

			while (!currentLevel.empty() || !neighborQueue.empty()) {
				if (currentLevel.empty()) {
					level = depth[neighborQueue.back()];
				}

				for (node u : currentLevel) {
					assert(depth[u] == level);
					processNode(u);
				}

				while (!neighborQueue.empty() && depth[neighborQueue.back()] == level) {
					node u = neighborQueue.back();
					neighborQueue.pop_back();
					assert(depth[u] == level);

					if (nodeTouched[u]) continue; // if the node was touched in the previous level, it was in currentLevel and thus has already been processed

					processNode(u);
				}

				--level;
				currentLevel.clear();
				currentLevel.swap(nextLevel);
			}

			count bestEdits = G.degree(nodeToMove) - rootMaxGain;

			if (rootEdits > rootMaxGain) {
				bestParent = none;
				bestEdits = G.degree(nodeToMove) - rootEdits;
			}

			std::vector<node> bestChildren;

			for (node u : touchedNodes) {
				if (u != nodeToMove && dynamicForest.parent(u) == bestParent && editDifference[u] != none && editDifference[u] > 0) {
					bestChildren.push_back(u);
				}
			}

			std::vector<count> missingBelow, missingAbove, existingBelow, existingAbove;

			bool countExactEdits = moveSubtrees;
#ifndef NDEBUG
			countExactEdits = true;
#endif

			if (countExactEdits) {
				missingBelow.resize(G.upperNodeIdBound(), 0);
				missingAbove.resize(G.upperNodeIdBound(), 0);
				existingBelow.resize(G.upperNodeIdBound(), 0);
				existingAbove.resize(G.upperNodeIdBound(), 0);

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

				assert(missingBelow[nodeToMove] == 0);
				assert(existingBelow[nodeToMove] == 0);

				for (node c : curChildren) {
					missingBelow[nodeToMove] += missingBelow[c];
					existingBelow[nodeToMove] += existingBelow[c];
				}

				if (curParent != none) {
					missingAbove[nodeToMove] = missingAbove[curParent];
					existingAbove[nodeToMove] = existingAbove[curParent];
				}
			}

#ifndef NDEBUG
			assert(curEdits == G.degree(nodeToMove) - existingAbove[nodeToMove] - existingBelow[nodeToMove] + missingAbove[nodeToMove] + missingBelow[nodeToMove]);

			count minEdits = curEdits;
			std::vector<node> minChildren;
			node minParent = curParent;

			G.forNodes([&](node u) {
				if (u == nodeToMove) return;
				if (existingBelow[u] >= missingBelow[u] || (editDifference[u] > 0 && editDifference[u] != none)) {
					assert(editDifference[u] == existingBelow[u] - missingBelow[u]);
				} else if (nodeTouched[u]) {
					assert(editDifference[u] == none);
				}
			});

			G.forNodes([&](node u) {
				if (dynamicForest.children(u).empty() && marker[u]) {
					assert(editDifference[u] == 1);
				}
			});

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
						if (dynamicForest.children(c).empty() && marker[c]) {
							assert(editDifference[c] == 1);
						}
						assert(editDifference[c] == existingBelow[c] - missingBelow[c]);

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

			dynamicForest.dfsFrom(none, [](node){}, tryEditBelow);
			tryEditBelow(none);

			assert(minEdits == bestEdits);

			count editDifferenceControl = G.degree(nodeToMove);
			if (bestParent != none) {
				editDifferenceControl -= (existingAbove[bestParent] - missingAbove[bestParent]);
			}

			for (node u : bestChildren) {
				editDifferenceControl -= (existingBelow[u] - missingBelow[u]);
			}
			assert(minEdits == editDifferenceControl);

			TRACE("Current edits: ", curEdits, " (with parent ", curParent, " and current children ", curChildren, "), minimum edits: ", minEdits);
			TRACE("Quadratic algorithm wants to have new parent ", minParent, " and new children ", minChildren);
			TRACE("Linear algorithm wants to have new parent ", bestParent, " and new children ", bestChildren);
#endif
			bool moveWithSubtree = false;

			// calculate the number of saved edits as comparing the absolute number of edits doesn't make sense
			count savedEdits = curEdits - bestEdits;

			// cleanup for linear move
			for (node u : touchedNodes) {
				lastVisitedDFSNode[u] = u;
				maxGain[u] = 0;
				editDifference[u] = 0;
				nodeTouched[u] = false;
				bestParentBelow[u] = none;
			}

			touchedNodes.clear();

			G.forEdgesOf(nodeToMove, [&](node v) {
				marker[v] = false;
			});

			if (moveSubtrees) {
				dynamicForest.setParent(nodeToMove, curParent);
				for (node c : curChildren) {
					dynamicForest.setParent(c, nodeToMove);
				}

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
				}, [](node) {});

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
					[](node) {});
				});

				// virtually remove the subtree of u from the forest concerning the missing/existingBelow-counters
				for (node u = curParent; u != none; u = dynamicForest.parent(u)) {
					existingBelow[u] -= existingBelow[nodeToMove];
					missingBelow[u] -= missingBelow[nodeToMove];
				}

				// calculate how many edits the whole subtree currently needs
				count curSubtreeEdits = subtreeExtDegree;
				if (curParent != none) { // FIXME here we ignore edits inside the tree as they do not change. Is this okay?
					curSubtreeEdits += missingAbove[curParent];
					curSubtreeEdits -= existingAbove[curParent];
				}

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
						bestEdits = edits;
						bestChildren = std::move(children);
						bestParent = p;
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

				TRACE("After subtree moving, ", savedEdits, " edits will be saved");
				TRACE("After subtree moving (quadratic) algorithm wants to have new parent ", bestParent, " and new children ", bestChildren);
			}


			if (savedEdits > 0) {
				if (!moveWithSubtree) {
					dynamicForest.isolate(nodeToMove);
				}
				dynamicForest.setParent(nodeToMove, bestParent);
				for (node c : bestChildren) {
					dynamicForest.setParent(c, nodeToMove);
				}

				hasMoved = true;
				numEdits -= savedEdits;

#ifndef NDEBUG
				forest = dynamicForest.toGraph();
				assert(numEdits == countNumberOfEdits());
#endif
			} else if (!moveSubtrees) {
				dynamicForest.setParent(nodeToMove, curParent);
				for (node c : curChildren) {
					dynamicForest.setParent(c, nodeToMove);
				}
			}

			dynamicForest.dfsFrom(nodeToMove, [&](node c) {
				if (dynamicForest.parent(c) != none) {
					depth[c] = depth[dynamicForest.parent(c)] + 1;
				} else {
					depth[c] = 0;
				}
			}, [](node){});
		});

		usedIterations = i+1;
	}

	forest = dynamicForest.toGraph();

	assert(numEdits == countNumberOfEdits());
}

NetworKit::count NetworKit::QuasiThresholdEditingLocalMover::countNumberOfEdits() const {
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

	return numMissingEdges + (G.numberOfEdges() - numExistingEdges);
}


std::vector< NetworKit::node > NetworKit::QuasiThresholdEditingLocalMover::getParents() const {
	std::vector<node> parents(G.upperNodeIdBound());

	G.forNodes([&](node u) {
		node p = none;
		forest.forNeighborsOf(u, [&](node v) {
			p = v;
		});

		parents[u] = p;
	});

	return parents;
}


NetworKit::count NetworKit::QuasiThresholdEditingLocalMover::getNumberOfEdits() const {
	return numEdits;
}

NetworKit::Graph NetworKit::QuasiThresholdEditingLocalMover::getQuasiThresholdGraph() const {
	TreeReachabilityGraphGenerator gen(forest);
	gen.run();
	return gen.getGraph();
}

NetworKit::count NetworKit::QuasiThresholdEditingLocalMover::getUsedIterations() const {
	return usedIterations;
}
