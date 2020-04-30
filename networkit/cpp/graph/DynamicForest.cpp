/*
 *
 */

#include <networkit/graph/DynamicForest.hpp>
#include <set>

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

		assert(nodes[u].parent == none);
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

void NetworKit::DynamicForest::moveUpNeighbor(node referenceNode, node neighbor) {
	if(path_membership[neighbor].valid()){
		TRACE("Move up ", neighbor);
		SimplePathPtr sp = path_membership[neighbor];
		if(sp->referenceNode != referenceNode){
			sp->referenceNode = referenceNode;
			sp->neighborCount = 0;
		}
		count oldPos = path_pos[neighbor];
		//neighbor already considered
		if(oldPos < sp->neighborCount){
			return;
		}
		node firstNonNeighbor = sp->pathNodes[sp->neighborCount];
		sp->pathNodes[sp->neighborCount] = neighbor;
		path_pos[neighbor] = sp->neighborCount;
		sp->pathNodes[oldPos] = firstNonNeighbor;
		path_pos[firstNonNeighbor] = oldPos;
		sp->neighborCount++;
		
		swapNodes(neighbor, firstNonNeighbor);
		
	}
	assert(pathsValid());
}

void NetworKit::DynamicForest::swapNodes(node u, node v){
	node oldParent = parent(u);
	std::vector<node> oldChildren = nodes[u].children;
	setParent(u, parent(v));
	for(node c : nodes[v].children){
		setParent(c, u);
	}
	setParent(v, oldParent);
	for(node c : oldChildren){
		setParent(c, v);
	}
	count oldDepth = d[u];
	d[u] = d[v];
	d[v] = oldDepth;
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
		if(oldPos == length - 1){
			return;
		}
		if(oldPos == length - 2){
			removeFromPath(oldPath->pathNodes[length -1]);
			return;
		}
		if (oldPos == 0){
			removeFromPath(oldPath->pathNodes[0]);
			return;
		}
		SimplePathPtr newPath(new SimplePath);
		newPath->neighborCount = 0;
		newPath->referenceNode = none;
		newPath->pathNodes = std::vector<node>();
		for(int i = 0; i <= oldPos; i++){
			node nodeToMove = oldPath->pathNodes[i];
			path_pos[nodeToMove] = i;
			path_membership[nodeToMove] = newPath;
			newPath->pathNodes.push_back(nodeToMove);
		}
		for(int i = oldPos + 1; i < length; i++){
			node nodeToMove = oldPath->pathNodes[i];
			oldPath->pathNodes[i - (oldPos + 1)] = nodeToMove;
			path_pos[nodeToMove] = i - (oldPos + 1);
		}
		for(int i = 0; i <= oldPos; i++){
			oldPath->pathNodes.pop_back();
		}
	}
	
	
}

void NetworKit::DynamicForest::unionPathWith(node moveNode, node keepNode){
	TRACE("Union paths of ", moveNode, " and ", keepNode);
	/*if(!path_membership[keepNode].valid()){
		SimplePathPtr newPath(new SimplePath);
		newPath->neighborCount = 0;
		newPath->referenceNode = none;
		newPath->pathNodes = std::vector<node>();
		newPath->pathNodes.push_back(keepNode);
		path_membership[keepNode] = newPath;
		path_pos[keepNode] = 0;
	}
	
	SimplePathPtr newPath = path_membership[keepNode];
	
	if(!path_membership[moveNode].valid()){
		path_pos[moveNode] = newPath->pathNodes.size();
		path_membership[moveNode] = newPath;
		newPath->pathNodes.push_back(moveNode);
		return;		
	}
	SimplePathPtr oldPath = path_membership[moveNode];
	count length = oldPath->pathNodes.size();
	for(int i = 0; i < length; i++){
		//hier müssen die Knoten eigentlich vorne 
	}*/
	if(moveNode == keepNode){
		return;
	}
	if(path_membership[moveNode].valid() && path_membership[keepNode].valid() && path_membership[moveNode] == path_membership[keepNode]){
		return;
	}
	SimplePathPtr newPath(new SimplePath);
	newPath->neighborCount = 0;
	newPath->referenceNode = none;
	newPath->pathNodes = std::vector<node>();
	newPath->pathNodes.size();
	count moveLength;
	if(!path_membership[moveNode].valid()){
		newPath->pathNodes.push_back(moveNode);
		path_membership[moveNode] = newPath;
		path_pos[moveNode] = 0;
		moveLength = 1;
	} else {
		SimplePathPtr movePath = path_membership[moveNode];
		moveLength = movePath->pathNodes.size();
		for(int j = 0; j < moveLength; j++){
			node curr = movePath->pathNodes[j];
			path_membership[curr] = newPath;
			path_pos[curr] = j;
			newPath->pathNodes.push_back(curr);
		}
	}
	if(!path_membership[keepNode].valid()){
		newPath->pathNodes.push_back(keepNode);
		path_membership[keepNode] = newPath;
		path_pos[keepNode] = moveLength;
	} else {
		SimplePathPtr keepPath = path_membership[keepNode];
		count keepLength = keepPath->pathNodes.size();
		for(int j = 0; j < keepLength; j++){
			node curr = keepPath->pathNodes[j];
			path_membership[curr] = newPath;
			path_pos[curr] = j + moveLength;
			newPath->pathNodes.push_back(curr);
		}
	}
	
}


void NetworKit::DynamicForest::moveToPosition(node u, node p, const std::vector<node> &children){
	//assert(allChildren adopted);
	node oldP = parent(u);
	std::vector<node> oldChildren = nodes[p].children;
	count oldChildCount = oldChildren.size();
	
	for(node c : children){
		assert(std::find(oldChildren.begin(), oldChildren.end(), c) != oldChildren.end());
	}
	
	setParent(u, p);
	for (node c : children) {
		setParent(c, u);
	}
	TRACE("Move ", u, " below parent ", p, " adopting children ", children);
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
	
	TRACE("New place");
	assert(pathsValid());
	

}


bool NetworKit::DynamicForest::pathsValid(){
	TRACE(printPaths());
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
	for (node u = 0; u < path_membership.size(); u++){
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