#ifndef QUASITHRESHOLDEDITINGLOCALMOVER_H
#define QUASITHRESHOLDEDITINGLOCALMOVER_H

#include <networkit/graph/Graph.hpp>
#include <networkit/graph/DynamicForest.hpp>
#include <networkit/structures/Cover.hpp>
#include <networkit/auxiliary/SignalHandling.hpp>
#include <networkit/auxiliary/Log.hpp>

namespace NetworKit {
	
	
	
	
	class BucketQueue { 
		count nextNode;
		count currentBucket;
		std::vector<node> nodes;
		//points to first element in the bucket
		std::vector<count> border;
		
		public :
		BucketQueue(count n = 0) : nodes(n), border(n), nextNode(none), currentBucket(none){
		}
		
		
		void fill(const std::vector<node> &elements, const DynamicForest &dynamicForest){
			count max_depth = 2 * elements.size();
			std::fill(border.begin(), border.begin() + std::min(max_depth + 1, border.size()), 0);
			for (node u : elements) {
				if(dynamicForest.depth(u) > max_depth) continue;
				border[dynamicForest.depth(u)] += 1;
			}
			for (int m = 1; m <= max_depth; m++){
				border[m] += border[m-1];
			}
			currentBucket = max_depth;
	    nextNode = none;
			for(int j = elements.size() - 1; j >= 0; j--){
				node u = elements[j];
				if(dynamicForest.depth(u) > max_depth){
					continue;
				}
				border[dynamicForest.depth(u)] -= 1;
	    	nodes[border[dynamicForest.depth(u)]] = u;
				nextNode += 1;
	  	}
		
			if(nextNode == none){
				return;
			} 

				
		}
		
		node next(){
			if(nextNode == none){
				return none;
			}
			node result = nodes[nextNode];
			while(nextNode < border[currentBucket]){
				currentBucket -= 1;
			}
			nextNode -= 1;
			return result;
		}
		
		void insertParent(node p){
			nextNode += 1;
			//first element of currentBucket
			count bucketBorder = border[currentBucket];
			node firstOfBucket = nodes[bucketBorder];
			nodes[nextNode] = firstOfBucket;
			nodes[bucketBorder] = p;
			border[currentBucket] += 1;
		}
		
		bool empty() {
			return nextNode == none;
		}
		
		std::string printQueue(){
			if(empty()){
				return "BucketQueue:";
			}
			std::stringstream ss;
			ss << "BucketQueue:";
			int j = 0;
			for (count i = 0; i <= nextNode; i++){
				while(border[j] == i){
					ss << "| ";
					j++;
				}
				ss << nodes[i] << " ";
			}
			return ss.str();
		}
		
		
	}; 

class QuasiThresholdEditingLocalMover {
public:
	QuasiThresholdEditingLocalMover(const NetworKit::Graph &G, 
		const std::vector< NetworKit::node > &parent = std::vector<NetworKit::node>(), 
		NetworKit::count maxIterations = 2,  
		bool sortPaths = true,
		bool randomness = false,
		const std::vector< NetworKit::node > &order = std::vector<NetworKit::node>());

	void run();

	Graph getQuasiThresholdGraph() const;
	count getNumberOfEdits() const;
	count getUsedIterations() const;
	count getPlateauSize() const;
	std::vector<node> getParents() const;
	Cover getCover(NetworKit::index mergeDepth) const;
private:
	count countNumberOfEdits() const;
	count editsIncidentTo(node u) const;
	void localMove(node nodeToMove);
	void processNode(node u, node nodeToMove);
	void compareWithQuadratic(node nodeToMove) const;
	bool randomBool(count options) const;
	
	const Graph& G;
	Graph forest;
	count maxIterations;
	count usedIterations;
	bool sortPaths;
	bool randomness;
	std::vector<NetworKit::node> order;
	bool insertRun;
	count numEdits;
	
	
	Aux::SignalHandler handler;
	DynamicForest dynamicForest;
	bool hasMoved;
	std::vector<bool> marker;

	std::vector<count> numNeighbors;
	NetworKit::BucketQueue bucketQueue;
	std::vector<node> neighbors;
	std::vector<node> currentLevel;
	std::vector<node> touchedNodes;
	std::vector<node> lastVisitedDFSNode;
	std::vector<node> bestParentBelow;

	std::vector<count> maxGain;
	std::vector<count> editDifference;
	std::vector<bool> nodeTouched;
	
	count bestParent;
	count rootMaxGain;
	count rootEdits;
	
	count bestEdits;
	count curEdits;
	node curParent;
	std::vector<node> curChildren;
	std::vector<node> bestChildren;
	
	std::vector<bool> existing;
	std::vector<count> equalBestParents;
	count rootEqualBestParents;
	
	count editsBefore;
	count plateauSize;
	
	
	
};

} // namespace NetworKit

#endif // QUASITHRESHOLDEDITINGLOCALMOVER_H
