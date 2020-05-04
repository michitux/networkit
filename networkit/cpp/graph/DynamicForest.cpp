/*
 *
 */

#include <networkit/graph/DynamicForest.hpp>
#include <set>

NetworKit::DynamicForest::DynamicForest(const NetworKit::Graph &G) : 
nodes(G.upperNodeIdBound()),
path_membership(G.upperNodeIdBound(), none),
path_pos(G.upperNodeIdBound(), 0),
d(G.upperNodeIdBound(), 0),
freeList(G.upperNodeIdBound()),
paths(G.upperNodeIdBound(), SimplePath()) {
	
	std::iota(path_membership.begin(), path_membership.end(), 0);
	G.forNodes([&](node u) {
		paths[u].pathNodes.push_back(u);
		
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
			//node has no siblings
			if(nodes[parent(c)].children.size() == 1){
				unionPathWith(parent(c), c);
			}
			
		}
	}, [](node){});
	assert(pathsValid());
	TRACE("Dynamic Forest constructed");
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

	if (oldP == p || u == p) return;
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
	removeFromPath(u);
	dfsFrom(u, [&](node c) {
		--d[c];
	}, [](node){});
	d[u] = 0;
	
	std::vector<node> oldChildren = nodes[u].children;
	node p = nodes[u].parent;

	setParent(u, none);

	for (node c : oldChildren) {
		setParent(c, p);
	}
	
	if(p != none && nodes[p].children.size() == 1){
		unionPathWith(p, nodes[p].children[0]);
	}

	
	TRACE("Isolate ", u);
	assert(pathsValid());
	
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

void NetworKit::DynamicForest::moveUpNeighbor(node neighbor, node referenceNode) {
	SimplePath* sp = path(neighbor);
	if(sp->length() > 1){
		TRACE("Move up ", neighbor);
		if(sp->referenceNode != referenceNode){
			sp->referenceNode = referenceNode;
			sp->neighborCount = 0;
		}
		if(sp->neighborCount >= sp->length()) return;
		count oldPos = path_pos[neighbor];
		count neighborPos = sp->length() - 1 - sp->neighborCount;
		if(oldPos <= neighborPos){
			sp->neighborCount++;
			if(oldPos ==  neighborPos) return;
			node firstNonNeighbor = sp->pathNodes[neighborPos];
			sp->pathNodes[neighborPos] = neighbor;
			path_pos[neighbor] = neighborPos;
			sp->pathNodes[oldPos] = firstNonNeighbor;
			path_pos[firstNonNeighbor] = oldPos;
			if(oldPos == 0){
				for(node c : nodes[neighbor].children){
					setParent(c, firstNonNeighbor);
				}
			} else {
				setParent(sp->pathNodes[oldPos - 1], firstNonNeighbor);
			} 
			
			if (neighborPos == sp->length() - 1){
				setParent(neighbor, parent(firstNonNeighbor));
			} else {
				setParent(neighbor, sp->pathNodes[neighborPos + 1]);
			}
			if(oldPos ==  neighborPos - 1){
				setParent(firstNonNeighbor, neighbor);
			} else {
				setParent(sp->pathNodes[neighborPos - 1], neighbor);
				setParent(firstNonNeighbor, sp->pathNodes[oldPos + 1]);	
			}
			count oldDepth = d[neighbor];
			d[neighbor] = d[firstNonNeighbor];
			d[firstNonNeighbor] = oldDepth;
		}
		assert(pathsValid());		
	}
}


void NetworKit::DynamicForest::removeFromPath(node u) {
	SimplePath* sp = path(u);
	count pos = path_pos[u];
	if(sp->length() > 1){
		addToPath(u, newPath());
		for(int i = pos + 1; i < sp->length(); i++){
			node v = sp->pathNodes[i];
			sp->pathNodes[i - 1] = v;
			path_pos[v] = i - 1;
		}
		sp->pathNodes.pop_back();
	}

}

void NetworKit::DynamicForest::addToPath(node u, index newId){
	TRACE("Adding ", u, " to path ", newId);
	SimplePath* sp = &(paths[newId]);
	path_membership[u] = newId;
	path_pos[u] = sp->length();
	sp->pathNodes.push_back(u);
}

//split after u
void NetworKit::DynamicForest::splitPath(node u){
	SimplePath* oldPath = path(u);
	if(oldPath->length() > 1){
		TRACE("Split path after ", u);
		count oldPos = path_pos[u];
		if (oldPos == 0){
			return;
		}
		if (oldPos == 1){
			removeFromPath(oldPath->pathNodes[0]);
			return;
		}
		if(oldPos == oldPath->length() - 1){
			removeFromPath(oldPath->pathNodes[oldPath->length() -1]);
			return;
		}
		index newId = newPath();
		SimplePath* newPath = &paths[newId];
		for(int i = oldPos; i < oldPath->length(); i++){
			addToPath(oldPath->pathNodes[i], newId);
		}
		oldPath->pathNodes.erase(oldPath->pathNodes.begin() + oldPos, oldPath->pathNodes.end());
	}
	
	
}

void NetworKit::DynamicForest::unionPathWith(node moveNode, node keepNode){
	TRACE("Union paths of ", moveNode, " and ", keepNode);
	if(path_membership[moveNode] == path_membership[keepNode]){
		return;
	}
	index keepPathId = path_membership[keepNode];
	index movePathId = path_membership[moveNode];
	SimplePath* movePath = path(moveNode);
	count l = movePath->length();
	for(count i = 0; i < l; i++){
		addToPath(movePath->pathNodes[i], keepPathId);
	}
	deletePath(movePathId);
}


void NetworKit::DynamicForest::moveToPosition(node u, node p, const std::vector<node> &children){
	node oldP = parent(u);
	count oldChildCount = none;
	if(p != none){
		std::vector<node> oldChildren = nodes[p].children;
		//TRACE("Move ", u, " below parent ", p, " adopting children ", children, " out of ", oldChildren);
		oldChildCount = oldChildren.size();
		for(node c : children){
			assert(std::find(oldChildren.begin(), oldChildren.end(), c) != oldChildren.end());
		}
	} else {
		//TRACE("Make ", u, " root adopting children ", children);
		for(node c : children){
			assert(std::find(roots.begin(), roots.end(), c) != roots.end());
		}
	}
	
	setParent(u, p);
	for (node c : children) {
		setParent(c, u);
	}
	removeFromPath(u);
	if(p != none){
		splitPath(p);
	}
	if(oldP != none && nodes[oldP].children.size() == 1){
		unionPathWith(oldP, nodes[oldP].children[1]);
	} 
	
	if(children.size() == 1){
		unionPathWith(u, children[0]);
	}
	if(p != none && oldChildCount == children.size()){
		unionPathWith(p, u);
	}
	
	dfsFrom(u, [&](node c) {
		if (parent(c) != none) {
			d[c] = d[parent(c)] + 1;
		} else {
			d[c] = 0;
		}
	}, [](node){});
	
	assert(pathsValid());
	

}


bool NetworKit::DynamicForest::pathsValid(){
	//TRACE(printPaths());
	//checking if tree is acyclic
	/*for(node u = 0; u < nodes.size(); u++){
		std::vector<bool> marked(nodes.size(), 0);
		dfsFrom(u, [&](node c) {
			if(marked[c] == 1){
				TRACE(c, " visited two times from ", u);
			}
			assert(marked[c] == 0);
			marked[c] = 1;
			}, [](node){});
	}*/

	for(int i = 0; i < nodes.size(); i++){
		TreeNode u = nodes[i];
		for(index j = 0; j < children(i).size(); j++){
			node c = children(i)[j];
			assert(parent(c) == i);
			assert(nodes[c].posInParent == j);
		}
	}
	
	dfsFrom(none, [&](node c) {
			if (c != none && parent(c) != none) {
				assert(d[c] == d[parent(c)] + 1);
				//node has no siblings
				if(nodes[parent(c)].children.size() == 1){
					assert(path_membership[c] == path_membership[parent(c)]);
				}
			}
		}, [](node){});
	for (node u = 0; u < nodes.size(); u++){
		std::vector<node> pNodes = path(u)->pathNodes;
		assert(pNodes[path_pos[u]] == u);
		for(int i = 0; i < pNodes.size(); i++){
			node v = pNodes[i];
			assert(path_membership[v] == path_membership[u]);
			assert(path_pos[v] == i);
			if(i > 0){
				assert(nodes[v].children.size() == 1 && nodes[v].children[0] == pNodes[i-1]);
			} else {
				assert(nodes[v].children.size() != 1);
			}
		}
	}
	return 1;
}


std::string NetworKit::DynamicForest::printPaths(){
	std::stringstream ss;
	for(node u = 0; u < path_membership.size(); u++){
		ss << "{" << u << " ";
		ss << "[";
		std::vector<node> pNodes = path(u)->pathNodes;
		for(node u : pNodes){
			ss << u << " ";
		}
		ss << "]";
		if(path_pos[u] != none){
			ss << " pos " << path_pos[u];
		}
		ss << "} ";
	}
	return ss.str();
}