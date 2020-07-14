#ifndef QUASITHRESHOLDEDITINGLOCALMOVER_H
#define QUASITHRESHOLDEDITINGLOCALMOVER_H

#include <networkit/community/QuasiThresholdMover/QuasiThresholdEditingLinear.hpp>
#include <networkit/community/QuasiThresholdMover/EditingRunner.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/graph/GraphTools.hpp>

namespace NetworKit {
	

class QuasiThresholdEditingLocalMover {
public:
	QuasiThresholdEditingLocalMover(const NetworKit::Graph &G, 
		Initialization initializarion = TRIVIAL, 
		NetworKit::count maxIterations = 5,  
		bool sortPaths = true,
		bool randomness = false,
		count maxPlateauSize = 4,
		bool useBucketQueue = true);

	void run();

	Graph getQuasiThresholdGraph() const;
	count getNumberOfEdits() const;
	count getUsedIterations() const;
	count getPlateauSize() const;
	count getRootEqualBestParents() const;
	
	void setInsertionOrder(std::vector<node> order);
private:

	const Graph& G;
	Initialization initialization;
	count maxIterations;
	bool sortPaths;
	bool randomness;
	count maxPlateauSize;
	bool insertRun;
	bool useBucketQueue;
	
	std::vector<node> order;
	
	count usedIterations;
	count numEdits;
	count plateauSize;
	count rootEqualBestParents;
	
	Graph quasiThresholdGraph;

};


} // namespace NetworKit

#endif // QUASITHRESHOLDEDITINGLOCALMOVER_H
