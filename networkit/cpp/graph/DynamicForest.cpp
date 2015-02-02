/*
 *
 */

#include "DynamicForest.h"

NetworKit::DynamicForest::DynamicForest(const NetworKit::Graph &G) : nodes(G.upperNodeIdBound()) {
	G.forNodes([&](node u) {
		if (G.degreeOut(u) == 0) {
			nodes[u].posInParent = roots.size();
			roots.push_back(u);
		}

		G.forEdgesOf(u, [&](node v) {
			nodes[u].parent = v;
			nodes[u].posInParent = nodes[v].children.size();
			nodes[v].children.push_back(u);
		});
	});
}

std::vector< NetworKit::node > NetworKit::DynamicForest::children(NetworKit::node u) const {
	std::vector<node> result;

	forChildrenOf(u, [&](node c) {
		result.emplace_back(c);
	});

	return result;
}


void NetworKit::DynamicForest::setParent(NetworKit::node u, NetworKit::node p) {
	node oldP = parent(u);

	if (oldP == p) return;

	index oldPos = nodes[u].posInParent;

	if (oldP != none) {
		nodes[oldP].children[oldPos] = nodes[oldP].children.back();
		nodes[oldP].children.pop_back();
		nodes[nodes[oldP].children[oldPos]].posInParent = oldPos;
	} else {
		roots[oldPos] = roots.back();
		roots.pop_back();
		nodes[roots[oldPos]].posInParent = oldPos;
	}

	nodes[u].parent = p;

	if (p == none) {
		nodes[u].posInParent = roots.size();
		roots.push_back(u);
	} else {
		nodes[u].posInParent = nodes[p].children.size();
		nodes[p].children.push_back(u);
	}
}

void NetworKit::DynamicForest::isolate(NetworKit::node u) {
	std::vector<node> oldChildren = nodes[u].children;
	node p = nodes[u].parent;

	setParent(u, none);

	for (node c : oldChildren) {
		setParent(c, p);
	}
}

NetworKit::Graph NetworKit::DynamicForest::toGraph() const {
	Graph result(nodes.size(), false, true);

	for (node r : roots) {
		dfsFrom(r, [](node u) {}, [&](node u) {
			if (parent(u) != none) {
				result.addEdge(u, parent(u));
			}
		});
	}

	return result;
}

