#ifndef DYNAMICFOREST_H
#define DYNAMICFOREST_H


#include <networkit/Globals.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {
	
	typedef long unsigned int pid;
	
	struct SimplePath {
		std::vector<node> pathNodes;
		count neighborCount;
		node referenceNode;
		pid parent;
		std::vector<pid> childPaths;
		index posInParent;
		count depth;
		
		SimplePath () : neighborCount(0) , referenceNode(none), parent(none), depth(0), posInParent(0){};
		
		count length() const {return pathNodes.size();};
		node upperEnd() const {return pathNodes.back();};
		node lowerEnd() const {return pathNodes[0];};
		
		void reset(){
			neighborCount = 0;
			referenceNode = none;
			parent = none;
			posInParent = 0;
			depth = 0;
			posInParent = 0;
		};
		
		std::string printPathInfo() const{
			std::stringstream ss;
			ss << "["<<"\n";
			ss << "pathNodes:";
			for(node x : pathNodes) ss << " " << x;
			ss << "\n";
			ss << "neighborCount: " << neighborCount << "\n";
			ss << "referenceNode: " << referenceNode << "\n";
			ss << "parent: " << parent << "\n";
			ss << "childPaths:";
			for(pid x : childPaths) ss << " " << x;
			ss << "\n";
			ss << "posInParent: " << posInParent << "\n";
			ss << "depth: " << depth << "\n";
			ss << "]" << "\n";			
			return ss.str();
		};

		
	};	
	

class DynamicForest {
public:
	DynamicForest();
	DynamicForest(const Graph& G);
	DynamicForest(std::vector<node> parents);
	node parent(node u) const;
	count depth(node u) const;
	std::vector<node> children(node u) const;
	Graph toGraph() const;
	void isolate(node u);
	void moveUpNeighbor(node referenceNode, node Neighbor);	
	void moveToPosition(node u, node p, const std::vector<node> &adoptedChildren);
	std::string printPaths() const;
	

	node nextChild(node child, node p) const {
		assert(parent(child) == p);
		if(isLowerEnd(p)){
			pid nextChildPath = nextPathChild(path(child), path(p));
			if(nextChildPath != path(p)) return paths[nextChildPath].upperEnd();
		} 
		return p;
	};
	
		

	node nextDFSNodeOnEnter(node curNode, node basis) const {
		if (curNode == none) {
			if (roots.empty()) return none; // forest is empty
			return paths[roots.front()].upperEnd();
		}

		if (childCount(curNode) > 0) {
			return children(curNode).front();
		} else if (curNode == basis) {
			return basis;
		} else {
			node u = curNode;
			node p = parent(u);

			assert(p == basis || p != none);

			while (p != basis && posInParent(u) == childCount(p) - 1) {
				u = p;
				p = parent(p);

				assert(p == basis || p != none);
			}

			index nextPos = posInParent(u) + 1;
			if (nextPos == childCount(p)) {
				return basis;
			} else {
				assert(isLowerEnd(p));
				return paths[paths[path(p)].childPaths[nextPos]].upperEnd();
			}
		}
	};
	
	
	template <typename F>
	void forChildrenOf(node u, F handle) const {
		if(u != none && !isLowerEnd(u)){
			handle(nextNodeInPath(u));
		} else {
			forPathChildrenOf(path(u), [&](pid c) {
				node child = paths[c].upperEnd();
				assert(child < path_membership.size());
				handle(paths[c].upperEnd());
			});
		}
	}



	template <typename F1, typename F2>
	void dfsFrom(node u, F1 onEnter, F2 onExit) const {
		struct DFSEvent {
			node n;
			bool isEnter;
			DFSEvent(node n, bool isEnter) : n(n), isEnter(isEnter) {};
		};

		std::stack<DFSEvent> toProcess;
		toProcess.emplace(u, false);
		toProcess.emplace(u, true);
		while (!toProcess.empty()) {
			DFSEvent ev = toProcess.top();
			toProcess.pop();

			if (ev.isEnter) {
				onEnter(ev.n);
				forChildrenOf(ev.n, [&](node c) {
					toProcess.emplace(c, false);
					toProcess.emplace(c, true);
				});
			} else {
				onExit(ev.n);
			}
		}
	}

	
private:
	
	pid path(node u) const;
	
	bool isUpperEnd (node u) const;
	bool isLowerEnd (node u) const;
	node nextNodeInPath(node u) const;
	node previousNodeInPath(node u) const;
	
	void deletePath(pid i);
	pid newPath();
	
	count childCount (node u) const;
	index posInParent (node u) const;

	void setParent(node u, node p);
	void setParentPath(pid s, pid p);
	void swapNodesWithinPath(node u, node v);
	void addToPath(node u, pid newId);
	void splitPath(pid sp, pid splitPos);
	void unionPaths(pid upperPath, pid lowerPath);
	void isolatePath(pid sp);
	void isolateNode(node u);
	void updateDepthInSubtree(pid start);
	
	bool pathsValid();
	
	std::vector<SimplePath> paths;
	std::vector<pid> freeList;
	std::vector<pid> path_membership;
	std::vector<index> path_pos;
	std::vector<pid> roots;
	
	
	
	pid nextPathChild(pid child, pid p) const{
		assert(paths[child].parent == p);

		index nextPos = paths[child].posInParent + 1;
		if (p == none) {
			if (nextPos < roots.size()) {
				return roots[nextPos];
			}
		} else {
			if (nextPos < paths[p].childPaths.size()) {
				return paths[p].childPaths[nextPos];
			}
		}
		return p;
	};
	
	
	template <typename F1, typename F2>
	void pathDfsFrom(pid u, F1 onEnter, F2 onExit) const {
		struct PathDFSEvent {
			pid sp;
			bool isEnter;
			PathDFSEvent(pid sp, bool isEnter) : sp(sp), isEnter(isEnter) {};
		};

		std::stack<PathDFSEvent> toProcess;
		toProcess.emplace(u, false);
		toProcess.emplace(u, true);

		while (!toProcess.empty()) {
			PathDFSEvent ev = toProcess.top();
			toProcess.pop();

			if (ev.isEnter) {
				onEnter(ev.sp);
				forPathChildrenOf(ev.sp, [&](pid c) {
					toProcess.emplace(c, false);
					toProcess.emplace(c, true);
				});
			} else {
				onExit(ev.sp);
			}
		}
	}

	template <typename F>
	void forPathChildrenOf(pid sp, F handle) const {
		if (sp == none) {
			for (pid r : roots) {
				handle(r);
			}
		} else {
			for (pid c : paths[sp].childPaths) {
				handle(c);
			}
		}
	}
	
};







} // namespace NetworKit

#endif // DYNAMICFOREST_H
