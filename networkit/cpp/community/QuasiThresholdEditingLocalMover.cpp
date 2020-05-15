/*
 *
 */

#include <unordered_set>

#include <networkit/community/QuasiThresholdEditingLocalMover.hpp>
#include <networkit/generators/TreeReachabilityGraphGenerator.hpp>
#include <networkit/graph/DynamicForest.hpp>
#include <networkit/graph/GraphTools.hpp>
#include <networkit/auxiliary/SignalHandling.hpp>
#include <networkit/auxiliary/Log.hpp>


NetworKit::QuasiThresholdEditingLocalMover::QuasiThresholdEditingLocalMover(
	const NetworKit::Graph &G, 
	const std::vector< NetworKit::node > &parent, 
	NetworKit::count maxIterations, 
	bool moveSubtrees,
	bool sortPaths,
	bool randomness,
	const std::vector<node> &order)
: G(G), 
	maxIterations(maxIterations), 
	moveSubtrees(moveSubtrees), 
	numEdits(none), 
	sortPaths(sortPaths),
	randomness(randomness),
	order(order) {
	forest = Graph(GraphTools::copyNodes(G), false, true);
	if(parent.size() == G.upperNodeIdBound()){
		G.forNodes([&](node u) {
			if (parent[u] != none) {
				forest.addEdge(u, parent[u]);
			}
		});
	}
	insertRun = (order.size() == G.upperNodeIdBound());
	dynamicForest = DynamicForest(forest);
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

	Aux::SignalHandler handler;
	handler.assureRunning();
	
 	

	numEdits = countNumberOfEdits();
	usedIterations = 0;

	hasMoved = true;
	marker = std::vector<bool>(G.upperNodeIdBound());

	

	bucketQueue = BucketQueue(G.upperNodeIdBound());
	lastVisitedDFSNode = std::vector<node>(G.upperNodeIdBound(), none);
	bestParentBelow = std::vector<node>(G.upperNodeIdBound(), none);
	G.parallelForNodes([&](node u) {
		lastVisitedDFSNode[u] = u;
	});
	maxGain = std::vector<count> (G.upperNodeIdBound(), 0);
	editDifference = std::vector<count> (G.upperNodeIdBound(), 0);
	nodeTouched  = std::vector<bool>(G.upperNodeIdBound(), false);
	
	equalBestParents = std::vector<count> (G.upperNodeIdBound(), 0);
	rootEqualBestParents = 0;
	
	count i;
	if(insertRun){
		existing = std::vector<bool>(G.upperNodeIdBound(), 0);
		i = 0;
	} else {
		existing = std::vector<bool>(G.upperNodeIdBound(), 1);
		i = 1;
	}

	for (; hasMoved && i <= maxIterations; ++i) {
		handler.assureRunning();
		hasMoved = false;

		if(insertRun){
			for(index j = 0; j < G.upperNodeIdBound();j++){
				node nodeToMove = order[j];
				localMove(nodeToMove);
				existing[nodeToMove] = 1;
			}
			insertRun = 0;
		} else {
			G.forNodesInRandomOrder([&](node nodeToMove) {
				localMove(nodeToMove);
			});
		}
		usedIterations = i;

		
	}

	forest = dynamicForest.toGraph();

	assert(numEdits == countNumberOfEdits());
}

void NetworKit::QuasiThresholdEditingLocalMover::localMove(node nodeToMove){
	handler.assureRunning();
	G.forEdgesOf(nodeToMove, [&](node v) {
		if(existing[v]){
			marker[v] = true;
			neighbors.push_back(v);
			if(sortPaths) {
				dynamicForest.moveUpNeighbor(v, nodeToMove);
			}
		}
	});
	if(sortPaths) {
		dynamicForest.moveUpNeighbor(nodeToMove, nodeToMove);
	}
	// remove the node from its tree but store the old position.
	curChildren = dynamicForest.children(nodeToMove);
	curParent = dynamicForest.parent(nodeToMove);
	curEdits = G.degree(nodeToMove);
	dynamicForest.dfsFrom(nodeToMove, [&](node c) {
		if (c != nodeToMove) {
			curEdits += 1 - 2 * marker[c];
		}
	}, [](node){});
	for (node p = dynamicForest.parent(nodeToMove); p != none; p = dynamicForest.parent(p)) {
		curEdits += 1 - 2 * marker[p];
	}
	dynamicForest.isolate(nodeToMove);
	bucketQueue.fill(neighbors, dynamicForest);
	bestParent = none;
	rootMaxGain = 0;
	rootEdits = 0;

	while(!bucketQueue.empty()){
		TRACE(bucketQueue.printQueue());
		node u = bucketQueue.next();
		processNode(u, nodeToMove);
	}

	bestEdits = G.degree(nodeToMove) - rootMaxGain;

	if (rootEdits > rootMaxGain) {
		bestParent = none;
		bestEdits = G.degree(nodeToMove) - rootEdits;
	}

	bestChildren.clear();
	for (node u : touchedNodes) {
		if (u != nodeToMove && dynamicForest.parent(u) == bestParent && editDifference[u] != none) {
			if(editDifference[u] > 0 || (randomness &&  randomBool(2))){
				bestChildren.push_back(u);
			}
		}
	}


#ifndef NDEBUG
	compareWithQuadratic(nodeToMove);
#endif

	// calculate the number of saved edits as comparing the absolute number of edits doesn't make sense
	count savedEdits = curEdits - bestEdits;

	// cleanup for linear move
	for (node u : touchedNodes) {
		equalBestParents[u] = 0;
		lastVisitedDFSNode[u] = u;
		maxGain[u] = 0;
		editDifference[u] = 0;
		nodeTouched[u] = false;
		bestParentBelow[u] = none;
	}

	neighbors.clear();
	touchedNodes.clear();
	rootEqualBestParents = 0;

	G.forEdgesOf(nodeToMove, [&](node v) {
		marker[v] = false;
	});

	if (savedEdits > 0 || (savedEdits == 0 && randomness && randomBool(2))) {
		dynamicForest.moveToPosition(nodeToMove, bestParent, bestChildren);

		hasMoved = true;
		numEdits -= savedEdits;

#ifndef NDEBUG
		forest = dynamicForest.toGraph();
		assert(numEdits == countNumberOfEdits());
#endif
	} else  {
		dynamicForest.moveToPosition(nodeToMove, curParent, curChildren);
	}
}


void NetworKit::QuasiThresholdEditingLocalMover::processNode(node u, node nodeToMove){
		TRACE("Processing node ", u, " of depth ", dynamicForest.depth(u), " (node to move: ", nodeToMove, ")");
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
		if (maxGain[u] > 0 
				|| (editDifference[u] != none && editDifference[u] != 0) 
				|| (randomness && editDifference[u] == 0  && randomBool(2))) {
			node p = dynamicForest.parent(u);

			if (p != none) {
				if (!nodeTouched[p]) {
					nodeTouched[p] = true;
					touchedNodes.push_back(p);
					if(!marker[p]){
						bucketQueue.insertParent(p);
					}
				}

				if (editDifference[u] != none) {
					assert(editDifference[u] <= maxGain[u]);
					editDifference[p] += editDifference[u];
				}
				bool coin = 0;
				if(randomness && (maxGain[u] == maxGain[p])){
					equalBestParents[p]++;
					coin = randomBool(equalBestParents[p]);
				}
				if (maxGain[u] > maxGain[p] || coin) {
					maxGain[p] = maxGain[u];
					bestParentBelow[p] = bestParentBelow[u];
					equalBestParents[p]=1;
				}
			} else {
				bool coin = 0;
				if(randomness && (maxGain[u] == rootMaxGain)){
					rootEqualBestParents++;
					coin = randomBool(rootEqualBestParents);
				}
				if (maxGain[u] > rootMaxGain || coin) {
					rootEqualBestParents = 1;
					rootMaxGain = maxGain[u];
					bestParent = bestParentBelow[u];
				}
				if (editDifference[u] != none) {
					rootEdits += editDifference[u];
				}
			}
		}

		if (dynamicForest.children(u).empty()) { assert(editDifference[u] == 1); }
}

void NetworKit::QuasiThresholdEditingLocalMover::compareWithQuadratic(node nodeToMove) const {
	std::vector<count> missingBelow, missingAbove, existingBelow, existingAbove;
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
			assert(curEdits == G.degree(nodeToMove) - existingAbove[nodeToMove] - existingBelow[nodeToMove] + missingAbove[nodeToMove] + missingBelow[nodeToMove]);

			count minEdits = curEdits;
			std::vector<node> minChildren;
			node minParent = curParent;
			G.forNodes([&](node u) {
				if (u == nodeToMove || dynamicForest.depth(u) > 2 * G.degree(nodeToMove)) return;
				if (existingBelow[u] >= missingBelow[u] || (editDifference[u] > 0 && editDifference[u] != none)) {
					assert(editDifference[u] == existingBelow[u] - missingBelow[u]);
				} else if (nodeTouched[u]) {
					assert(editDifference[u] == none);
				}
			});

			G.forNodes([&](node u) {
				if (dynamicForest.children(u).empty() && marker[u] && dynamicForest.depth(u) < 2*G.degree(nodeToMove)) {
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
					if (existingBelow[c] > missingBelow[c] && dynamicForest.depth(c) < 2*G.degree(nodeToMove)) { // TODO try >= (more children...)
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
}


NetworKit::count NetworKit::QuasiThresholdEditingLocalMover::countNumberOfEdits() const {
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


NetworKit::Cover NetworKit::QuasiThresholdEditingLocalMover::getCover(index mergeDepth) const {
	Cover c(G.upperNodeIdBound());

	DynamicForest dynamicForest(forest);

	std::vector<count> depth(G.upperNodeIdBound());

	index curSubset = none;

	dynamicForest.dfsFrom(none, [&](node u) {
		if (u != none) {
			if (dynamicForest.parent(u) == none) {
				depth[u] = 0;
			} else {
				depth[u] = depth[dynamicForest.parent(u)] + 1;
			}

			if (depth[u] == mergeDepth || (depth[u] < mergeDepth && dynamicForest.nextDFSNodeOnEnter(u, u) == u)) {
				curSubset = c.toSingleton(u);

				node p = dynamicForest.parent(u);
				while (p != none) {
					c.addToSubset(curSubset, p);
					p = dynamicForest.parent(p);
				}
			} else if (depth[u] > mergeDepth) {
				c.addToSubset(curSubset, u);
			}
		}
	}, [&](node) {
	});

	return c;
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

bool NetworKit::QuasiThresholdEditingLocalMover::randomBool(count options) const {
	return rand() % options == 0;
}
