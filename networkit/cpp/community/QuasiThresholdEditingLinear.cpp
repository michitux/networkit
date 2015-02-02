/*
 *
 */

#include "QuasiThresholdEditingLinear.h"
#include "../backbones/ChibaNishizekiTriangleCounter.h"

NetworKit::QuasiThresholdEditingLinear::QuasiThresholdEditingLinear(const NetworKit::Graph &G) : G(G), hasRun(false) {
}

void NetworKit::QuasiThresholdEditingLinear::run() {
	parent.clear();
	parent.resize(G.upperNodeIdBound(), none);

	std::vector<count> numConnectedAnchestors(G.upperNodeIdBound(), 0), parentStrength(G.upperNodeIdBound(), 0);

	ChibaNishizekiTriangleCounter triangleCounter;
	std::vector<int> triangles = triangleCounter.getAttribute(G, std::vector<int>());

	std::vector<count> pseudoP4C4(G.upperEdgeIdBound(), 0);

	G.parallelForEdges([&](node u, node v, edgeid eid) {
		pseudoP4C4[eid] = (G.degree(u) - 1 - triangles[eid]) * (G.degree(v) - 1 - triangles[eid]);
	});

	// bucket sort
	count n = G.numberOfNodes();
	std::vector<node> sortedNodes(n);
	{
		std::vector<index> nodePos(n + 1, 0);

		G.forNodes([&](node u) {
			++nodePos[n - G.degree(u)];
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
			sortedNodes[nodePos[n - G.degree(u)]++] = u;
		});
	}

	std::vector<bool> processed(G.upperNodeIdBound(), false);

	for (node u : sortedNodes) {
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
			}
		});
	}

	hasRun = true;

}


std::vector< NetworKit::node > NetworKit::QuasiThresholdEditingLinear::getParents() const {
	return parent;
}


NetworKit::Graph NetworKit::QuasiThresholdEditingLinear::getDefiningForest() const {
	Graph result(G.copyNodes(), false, true);

	if (!hasRun) throw std::runtime_error("Error, run must be called first");

	G.forNodes([&](node u) {
		if (parent[u] != none) {
			result.addEdge(u, parent[u]);
		}
	});

	return result;
}


NetworKit::Graph NetworKit::QuasiThresholdEditingLinear::getQuasiThresholdGraph() const {
	Graph result(G.copyNodes());

	if (!hasRun) throw std::runtime_error("Error, run must be called first");

	G.forNodes([&](node u) {
		for (node p = parent[u]; p != none; p = parent[p]) {
			result.addEdge(u, p);
		}
	});

	return result;
}

