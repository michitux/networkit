#ifndef DYNAMICFOREST_H
#define DYNAMICFOREST_H

#include "../Globals.h"
#include "Graph.h"
#include <unordered_set>
#include <stack>

namespace NetworKit {

class DynamicForest {
public:
	DynamicForest(const Graph& G);

	node parent(node u) const { return nodes[u].parent; };
	std::vector<node> children(node u) const;
	Graph toGraph() const;

	void setParent(node u, node p);
	void isolate(node u);

	template <typename F1, typename F2>
	void dfsFrom(node u, F1 onEnter, F2 onExit) const;

	node nextChild(node child, node parent) const {
		assert(nodes[child].parent == parent);

		index nextPos = nodes[child].posInParent + 1;
		if (parent == none) {
			if (nextPos < roots.size()) {
				return roots[nextPos];
			}
		} else {
			if (nextPos < nodes[parent].children.size()) {
				return nodes[parent].children[nextPos];
			}
		}

		return parent;
	};

	node nextDFSNodeOnEnter(node curNode, node basis) const {
		if (curNode == none) {
			if (roots.empty()) return none; // forest is empty
			return roots.front();
		}

		if (!nodes[curNode].children.empty()) {
			return nodes[curNode].children.front();
		} else if (curNode == basis) {
			return basis;
		} else {
			node u = curNode;
			node p = parent(u);

			assert(p == basis || p != none);

			while (p != basis && nodes[u].posInParent == nodes[p].children.size() - 1) {
				u = p;
				p = parent(p);

				assert(p == basis || p != none);
			}

			index nextPos = nodes[u].posInParent + 1;
			if (nextPos == nodes[p].children.size()) {
				return basis;
			} else {
				return nodes[p].children[nextPos];
			}
		}
	};

	template <typename F>
	void forChildrenOf(node u, F handle) const;
private:
	struct TreeNode {
		node parent;
		index posInParent;
		std::vector<node> children;
		TreeNode() : parent(none) {};
	};

	std::vector<node> roots;
	std::vector<TreeNode> nodes;
};

template <typename F1, typename F2>
void DynamicForest::dfsFrom(node u, F1 onEnter, F2 onExit) const {
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

template <typename F>
void DynamicForest::forChildrenOf(node u, F handle) const {
	if (u == none) {
		for (node r : roots) {
			handle(r);
		}
	} else {
		for (node c : nodes[u].children) {
			handle(c);
		}
	}
}



} // namespace NetworKit

#endif // DYNAMICFOREST_H
