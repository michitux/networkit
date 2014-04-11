/*
 * Dijkstra.h
 *
 *  Created on: Jul 23, 2013
 *      Author: Henning, Christian Staudt
 */

#ifndef DIJKSTRA_H_
#define DIJKSTRA_H_

#include "Graph.h"
#include "../auxiliary/PrioQueue.h"

namespace NetworKit {

/**
 * Dijkstra's SSSP algorithm.
 */
class Dijkstra {
protected:


public:

	Dijkstra(const Graph& G, node source);
	virtual ~Dijkstra();


	virtual void run();

	/** 
	 * return Vector of weighted distances from node @a source, i.e. the
 	 * length of the shortest path from @a source to any other vertex.
	 */
	virtual std::vector<edgeweight> getDistances() const;

	/**
	 * @return a shortest path from source node to target node @a t.
	 * Returns empty path if source and target are not connected.
	 */
	virtual std::vector<node> getPath(node t, bool forward=true) const;

	/**
	 * @return all shortest paths from source node to target node @a t.
	 * Returns empty set if source and target are not connected.
	 */
	virtual std::set<std::vector<node> > getPaths(node t, bool forward=true) const;

private:

	const Graph& G;
	const node source;
	std::vector<edgeweight> distances;
	std::vector<std::vector<node> > previous; // predecessors on shortest path
};

} /* namespace NetworKit */
#endif /* DIJKSTRA_H_ */
