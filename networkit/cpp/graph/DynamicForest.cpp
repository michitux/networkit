/*
 *
 */

#include <networkit/graph/DynamicForest.hpp>
#include <set>

NetworKit::DynamicForest::DynamicForest(const NetworKit::Graph &G) : 
nodes(G.upperNodeIdBound()),
path_membership(G.upperNodeIdBound(), nullptr),
path_pos(G.upperNodeIdBound(), none),
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
	if(path_membership[neighbor].valid()){
		TRACE("Move up ", neighbor);
		SimplePathPtr sp = path_membership[neighbor];
		if(sp->referenceNode != referenceNode){
			sp->referenceNode = referenceNode;
			sp->neighborCount = 0;
		}
		if(sp->neighborCount >= sp->pathNodes.size()) return;
		count oldPos = path_pos[neighbor];
		count neighborPos = sp->pathNodes.size() - 1 - sp->neighborCount;
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
			
			if (neighborPos == sp->pathNodes.size() - 1){
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
	if(path_membership[u].valid()){
		SimplePathPtr sp = path_membership[u];
		count pos = path_pos[u];
		path_membership[u] = nullptr;
		path_pos[u] = none;
		count length = sp->pathNodes.size();
		if(length == 2){
			node other;
			if(pos == 0){
				other = sp->pathNodes[1];
			} else {
				other = sp->pathNodes[0];
			}
			path_membership[other] = nullptr;
			path_pos[other] = none;
		} else {
			for(int i = pos + 1; i < length; i++){
				node v = sp->pathNodes[i];
				sp->pathNodes[i - 1] = v;
				path_pos[v] = i - 1;
			}
			sp->pathNodes.pop_back();
		}

	}

}

//split after u
void NetworKit::DynamicForest::splitPath(node u){
	if(path_membership[u].valid()){
		SimplePathPtr oldPath = path_membership[u];
		TRACE("Split path ", oldPath ," after ", u);
		count oldPos = path_pos[u];
		count length = oldPath->pathNodes.size();
		if (oldPos == 0){
			return;
		}
		if (oldPos == 1){
			removeFromPath(oldPath->pathNodes[0]);
			return;
		}
		if(oldPos == length - 1){
			removeFromPath(oldPath->pathNodes[length -1]);
			return;
		}
		SimplePathPtr newPath(new SimplePath);
		newPath->pathNodes = std::vector<node>(length - oldPos);
		for(int i = length - oldPos - 1; i >= 0 ; i--){
			node nodeToMove = oldPath->pathNodes[i + oldPos];
			newPath->pathNodes[i] = nodeToMove;
			path_pos[nodeToMove] = i;
			path_membership[nodeToMove] = newPath;
			oldPath->pathNodes.pop_back();
		}
	}
	
	
}

void NetworKit::DynamicForest::unionPathWith(node moveNode, node keepNode){
	TRACE("Union paths of ", moveNode, " and ", keepNode);
	if(moveNode == keepNode){
		return;
	}
	if(path_membership[moveNode].valid() && path_membership[keepNode].valid() && path_membership[moveNode] == path_membership[keepNode]){
		return;
	}
	
	if(!path_membership[keepNode].valid()){
		SimplePathPtr newPath(new SimplePath);
		newPath->pathNodes.push_back(keepNode);
		path_membership[keepNode] = newPath;
		path_pos[keepNode] = 0;
	}
	SimplePathPtr keepPath = path_membership[keepNode];
	count offset = keepPath->pathNodes.size();
	if(!path_membership[moveNode].valid()){
		keepPath->pathNodes.push_back(moveNode);
		path_membership[moveNode] = keepPath;
		path_pos[moveNode] = offset;
	} else {
		SimplePathPtr movePath = path_membership[moveNode];
		count moveLength = movePath->pathNodes.size();
		for(count i = 0; i < moveLength; i++){
			node currentNode = movePath->pathNodes[i];
			keepPath->pathNodes.push_back(currentNode);
			path_membership[currentNode] = keepPath;
			path_pos[currentNode] = offset;
			offset++;
		}

}
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
		for(node c : u.children){
			assert(parent(c) == i);
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
		if(!path_membership[u].valid()){
			assert(path_pos[u] == none);
		} else {
			std::vector<node> pNodes = path_membership[u]->pathNodes;
			assert(pNodes[path_pos[u]] == u);
			assert(pNodes.size() >= 2);
			for(int i = 0; i < pNodes.size(); i++){
				node v = pNodes[i];
				assert(path_membership[v] == path_membership[u]);
				assert(path_pos[v] == i);
				if(i > 0){
					assert(nodes[v].children.size() == 1 && nodes[v].children[0] == pNodes[i-1]);
				}
			}
		}
	}
	return 1;
}


std::string NetworKit::DynamicForest::printPaths(){
	std::stringstream ss;
	for(node u = 0; u < path_membership.size(); u++){
		ss << "{" << u << " ";
		if(path_membership[u].valid()){
			ss << "[";
			std::vector<node> pNodes = path_membership[u]->pathNodes;
			for(node u : pNodes){
				ss << u << " ";
			}
			ss << "]";
		} 
		if(path_pos[u] != none){
			ss << " pos " << path_pos[u];
		}
		ss << "} ";
	}
	return ss.str();
}