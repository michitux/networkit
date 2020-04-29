/*
 *
 */

#include <networkit/graph/DynamicForest.hpp>

NetworKit::DynamicForest::DynamicForest(const NetworKit::Graph &G) : 
nodes(G.upperNodeIdBound()),
path_membership(G.upperNodeIdBound(), nullptr),
path_pos(G.upperNodeIdBound(), none),
simplePaths(),
d(G.upperNodeIdBound(), 0) {
	
	G.forNodes([&](node u) {
		if (G.degreeOut(u) == 0) {
			nodes[u].posInParent = roots.size();
			roots.push_back(u);
			d[u] = 0;
		}

		G.forEdgesOf(u, [&](node v) {
			nodes[u].parent = v;
			nodes[u].posInParent = nodes[v].children.size();
			nodes[v].children.push_back(u);
		});
	});
	
dfsFrom(none, [&](node c) {
		if (c != none && parent(c) != none) {
			d[c] = d[parent(c)] + 1;
		}
	}, [](node){});
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
		/*if(nodes[oldP].children.size() == 1){
			unionPathWith(oldP, nodes[oldP].children[0]);
		}*/
	} else {
		roots[oldPos] = roots.back();
		roots.pop_back();
		nodes[roots[oldPos]].posInParent = oldPos;
	}

	nodes[u].parent = p;
	//splitPath(p);

	if (p == none) {
		nodes[u].posInParent = roots.size();
		roots.push_back(u);
	} else {
		nodes[u].posInParent = nodes[p].children.size();
		nodes[p].children.push_back(u);
		/*if(nodes[p].children.size() == 1){
			unionPathWith(p, nodes[p].children[0]);
		}*/
	}
}

void NetworKit::DynamicForest::isolate(NetworKit::node u) {
	//removeFromPath(u);
	std::vector<node> oldChildren = nodes[u].children;
	node p = nodes[u].parent;

	setParent(u, none);

	for (node c : oldChildren) {
		setParent(c, p);
	}
	
	d[u] = 0;
	dfsFrom(u, [&](node c) {
		--d[c];
	}, [](node){});
	
}

NetworKit::Graph NetworKit::DynamicForest::toGraph() const {
	Graph result(nodes.size(), false, true);

	for (node r : roots) {
		dfsFrom(r, [](node) {}, [&](node u) {
			if (parent(u) != none) {
				result.addEdge(u, parent(u));
			}
		});
	}

	return result;
}

void NetworKit::DynamicForest::moveUpNeighbor(node referenceNode, node Neighbor) {
	SimplePathPtr sp;
}

void NetworKit::DynamicForest::removeFromPath(node u) {

}

void NetworKit::DynamicForest::splitPath(node u){
	
}

void NetworKit::DynamicForest::unionPathWith(node moveNode, node keepNode){
	
}

void NetworKit::DynamicForest::moveToPosition(node u, node p, const std::vector<node> &children){
	//assert(allChildren adopted);
	node oldP = parent(u);
	count oldChildCount = nodes[p].children.size();
	
	setParent(u, p);
	for (node c : children) {
		setParent(c, u);
	}
	
	/*removeFromPath(u);
	splitPath(p);
	
	if(nodes[oldP].children.size() == 1){
		unionPathWith(oldP, nodes[oldP].children[1]);
	} 
	
	if(children.size() == 1){
		unionPathWith(children[0], u);
	}
	if(p != none && oldChildCount == children.size()){
		unionPathWith(p, u);
	}*/
	
	dfsFrom(u, [&](node c) {
		if (parent(c) != none) {
			d[c] = d[parent(c)] + 1;
		} else {
			d[c] = 0;
		}
	}, [](node){});
	

}


