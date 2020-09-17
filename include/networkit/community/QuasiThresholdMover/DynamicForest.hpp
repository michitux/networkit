#ifndef DYNAMICFOREST_H
#define DYNAMICFOREST_H


#include <networkit/Globals.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {

namespace QuasiThresholdMoving {
	class DynamicForest {
		
		
	public:
		DynamicForest();
		DynamicForest(const Graph& G, const std::vector<node>& parents = std::vector<node>());
		DynamicForest(const std::vector<node>& parents);
		
		void isolate(node u);
		void moveUpNeighbor(node referenceNode, node Neighbor);	
		void moveToPosition(node u, node p, const std::vector<node> &adoptedChildren);
		Graph toGraph() const;
		std::string printPaths() const;
		
		node parent(node u) const {
			if(u == none){
				return none;
			} else if(isUpperEnd(u)){
				pid p = paths[path(u)].parent;
				return (p == none) ? none : paths[p].lowerEnd();
			} else {
				return previousNodeInPath(u);
			}
		};

		/**
		 * Call the given callback for all ancestors of u excluding the virtual root.
		 *
		 * @param u        The node for whose ancestors @a callback shall be called.
		 * @param callback The callback to call
		 */
		template <typename F>
		void forAncestors(node u, F callback) const {
			assert(u != none);
			pid cp = path(u);
			index pp = path_pos[u] + 1;
			while (cp != none) {
				const auto& pn = paths[cp].pathNodes;
				for (auto it = pn.begin() + pp; it != pn.end(); ++it) {
					callback(*it);
				}

				cp = paths[cp].parent;
				pp = 0;
			}
		}

		count depth(node u) const {
			if(u ==  none){
				return none;
			} else {
				return paths[path(u)].depth + paths[path(u)].length() - 1 - path_pos[u];
			}
		};
		
		std::vector<node> children(node u) const {
			std::vector<node> result;
			forChildrenOf(u, [&](node c) {
				result.emplace_back(c);
			});
			
			return result;
		};
		
		
		
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
			if(!isLowerEnd(curNode)){
				return nextNodeInPath(curNode);
			} 
			if(!paths[path(curNode)].childPaths.empty()){
				return paths[paths[path(curNode)].childPaths.front()].upperEnd();
			} 
			if (curNode == basis) {
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
					handle(child);
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
		using pid = index;
		
		class SimplePath {
		public:
			SimplePath () : neighborCount(0) , referenceNode(none), parent(none), depth(0), posInParent(0){};
			count length() const {return pathNodes.size();};
			node upperEnd() const {return pathNodes.empty() ? none : pathNodes.back();};
			node lowerEnd() const {return pathNodes.empty() ? none : pathNodes[0];};
			void reset(){
				neighborCount = 0;
				referenceNode = none;
				parent = none;
				posInParent = 0;
				depth = 0;
				posInParent = 0;
			};
			
			
			pid parent;
			std::vector<pid> childPaths;
			index posInParent;
			std::vector<node> pathNodes;
			count neighborCount;
			node referenceNode;
			count depth;
			
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

		pid path(node u) const {
			if(u == none){
				return none;
			} else {
				assert(u < path_membership.size());
				return path_membership[u];
			}
		};

		bool isUpperEnd(node u) const {
			return (path_pos[u] == paths[path(u)].length() - 1);
		};
		bool isLowerEnd(node u) const {
			return (path_pos[u] == 0);
		};
		node nextNodeInPath(node u) const {
			if(isLowerEnd(u)){
				return none;
			} else {
				return paths[path(u)].pathNodes[path_pos[u] - 1];
			}
		};
		node previousNodeInPath(node u) const {
			if(isUpperEnd(u)){
				return none;
			} else {
				return paths[path(u)].pathNodes[path_pos[u] + 1];
			}
		};
		
		void addToPath(node u, pid newId){
			path_membership[u] = newId;
			path_pos[u] = paths[newId].length();
			paths[newId].pathNodes.push_back(u);
		};
		
		void setParentPath(pid s, pid p){
			pid oldP = paths[s].parent;
			if (oldP == p || s == p) return;
			index oldPos = paths[s].posInParent;
			//tidy up at old position
			if (oldP != none) {
				paths[oldP].childPaths[oldPos] = paths[oldP].childPaths.back();
				paths[oldP].childPaths.pop_back();
				paths[paths[oldP].childPaths[oldPos]].posInParent = oldPos;
			} else {
				roots[oldPos] = roots.back();
				roots.pop_back();
				paths[roots[oldPos]].posInParent = oldPos;
			}
			
			//place at new position
			paths[s].parent = p;
			if (p == none) {
				paths[s].posInParent = roots.size();
				roots.push_back(s);
			} else {
				paths[s].posInParent = paths[p].childPaths.size();
				paths[p].childPaths.push_back(s);
			}
		};

		void deletePath(pid i);
		pid newPath();
		
		count childCount (node u) const;
		index posInParent (node u) const;
		
		void setParent(node u, node p);
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
	
	
	
	
	
} //namespace QuasiThresholdMoving



} // namespace NetworKit

#endif // DYNAMICFOREST_H
