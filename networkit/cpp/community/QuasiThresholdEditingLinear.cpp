/*
 *
 */

#include <unordered_map>

#include <networkit/community/QuasiThresholdEditingLinear.hpp>
#include <networkit/edgescores/TriangleEdgeScore.hpp>
#include <networkit/auxiliary/SignalHandling.hpp>
#include <networkit/graph/GraphTools.hpp>

NetworKit::QuasiThresholdEditingLinear::QuasiThresholdEditingLinear(const NetworKit::Graph &G) : G(G), hasRun(false) {
}

void NetworKit::QuasiThresholdEditingLinear::run() {
	Aux::SignalHandler handler;
	parent.clear();
	parent.resize(G.upperNodeIdBound(), none);

	std::vector<count> numConnectedAnchestors(G.upperNodeIdBound(), 0), parentStrength(G.upperNodeIdBound(), 0);

	TriangleEdgeScore triangleCounter(G);
	triangleCounter.run();
	std::vector<count> triangles = triangleCounter.scores();

	handler.assureRunning();

	std::vector<count> pseudoP4C4(G.upperEdgeIdBound(), 0);

	G.parallelForEdges([&](node u, node v, edgeid eid) {
		pseudoP4C4[eid] = (G.degree(u) - 1 - triangles[eid]) * (G.degree(v) - 1 - triangles[eid]);
	});

	// bucket sort
	count n = G.numberOfNodes();
	std::vector<node> sortedNodes(n);
	std::vector<count> degree(G.upperNodeIdBound());
	std::vector<index> sortedPos(G.upperNodeIdBound());
	std::vector<index> nodePos(n + 1, 0);

	G.forNodes([&](node u) {
		degree[u] = G.degree(u);
		++nodePos[n - degree[u]];
	});

	// exclusive prefix sum
	index tmp = nodePos[0];
	index sum = tmp;
	nodePos[0] = 0;

	for (index i = 1; i < nodePos.size(); ++i) {
		tmp = nodePos[i];
		nodePos[i] = sum;
		sum += tmp;
	}

	G.forNodes([&](node u) {
		sortedPos[u] = nodePos[n - degree[u]]++;
		sortedNodes[sortedPos[u]] = u;
	});

	handler.assureRunning();

	std::vector<bool> processed(G.upperNodeIdBound(), false);

	for (node u : sortedNodes) {
		handler.assureRunning();
		processed[u] = true;

		std::unordered_map<node, count> parents;

		G.forNeighborsOf(u, [&](node, node v, edgeid eid) {
			// node v has already been removed
			if (processed[v]) return;

			parents[parent[v]] += (parent[v] == parent[u] || (pseudoP4C4[eid] <= parentStrength[v] && numConnectedAnchestors[v] <= (count)triangles[eid] + 1));
		});

		count maxParent = parent[u], ownParent = parent[u];
		count maxParentCount = parents[maxParent];

		for (auto parentCount : parents) {
			if (parentCount.second > maxParentCount) {
				maxParent = parentCount.first;
				maxParentCount = parentCount.second;
			}
		}

		if (maxParent != ownParent) {
			parent[u] = maxParent;
			parentStrength[u] = none;
			numConnectedAnchestors[u] = 0;
		}

		G.forNeighborsOf(u, [&](node, node v, edgeid eid) {
			if (processed[v]) return;

			if (parent[v] == parent[u] || (
				pseudoP4C4[eid] < parentStrength[v] && numConnectedAnchestors[v] < (count)triangles[eid] + 1
			)) {
				parent[v] = u;
				parentStrength[v] = pseudoP4C4[eid];
				numConnectedAnchestors[v] += 1;
			} else if (degree[v] > 0) { // we have ignored an edge of v - decrease degree of v by 1
				index newPos = nodePos[n - degree[v]] - 1; // this is the last entry of the bucket where v is
				index oldPos = sortedPos[v];
				assert(newPos >= oldPos);
				node x = sortedNodes[newPos];
				assert(degree[x] == degree[v]);
				std::swap(sortedNodes[oldPos], sortedNodes[newPos]);
				--nodePos[n - degree[v]];
				--degree[v];
				if (nodePos[n - degree[v]] < n) {
					assert(degree[sortedNodes[nodePos[n - degree[v]]]] < degree[v]);
				}
				if (nodePos[n - degree[x]] < n) {
					assert(degree[sortedNodes[nodePos[n - degree[x]]]] < degree[x]);
				}
				std::swap(sortedPos[v], sortedPos[x]);
				#ifndef NDEBUG
				count oldDeg = std::numeric_limits<count>::max();
				for (node y : sortedNodes) {
					assert(oldDeg >= degree[y]);
					oldDeg = degree[y];
				}
				#endif
			}
		});
	}

	hasRun = true;

}


std::vector< NetworKit::node > NetworKit::QuasiThresholdEditingLinear::getParents() const {
	return parent;
}


NetworKit::Graph NetworKit::QuasiThresholdEditingLinear::getDefiningForest() const {
	Graph result(GraphTools::copyNodes(G), false, true);

	if (!hasRun) throw std::runtime_error("Error, run must be called first");

	G.forNodes([&](node u) {
		if (parent[u] != none) {
			result.addEdge(u, parent[u]);
		}
	});

	return result;
}


NetworKit::Graph NetworKit::QuasiThresholdEditingLinear::getQuasiThresholdGraph() const {
	Graph result(GraphTools::copyNodes(G));

	if (!hasRun) throw std::runtime_error("Error, run must be called first");

	G.forNodes([&](node u) {
		for (node p = parent[u]; p != none; p = parent[p]) {
			result.addEdge(u, p);
		}
	});

	return result;
}

