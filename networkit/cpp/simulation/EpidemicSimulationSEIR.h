/*
 * EpidemicSimulationSEIR.h
 *
 *  Created on: 20.11.2015
 *      Author: Christian Staudt
 */

#ifndef EPIDEMICSIMULATIONSEIR_H_
#define EPIDEMICSIMULATIONSEIR_H_

#include "../base/Algorithm.h"
#include "../graph/Graph.h"

namespace NetworKit {

/**
 * @ingroup simulation
 * Node centrality index which ranks nodes by their degree.
 * Optional normalization by maximum degree.
 */
class EpidemicSimulationSEIR: public NetworKit::Algorithm {
public:
	/**
	 *
	 *
	 * @param G The network.
	 *
	 */
	EpidemicSimulationSEIR(const Graph& G, count tMax, double transP, count eTime, count iTime);

	void run() override;

protected:

	const Graph& G;
	count tMax;
	double transP;
	count eTime;
	count iTime;
	enum class State {S, E, I, R, U}; // Susceptible, Exposed, Infectious, Removed, Undefined
	std::vector<State> state;
	std::vector<index> timestamp;

};

} /* namespace NetworKit */

#endif /* DEPIDEMICSIMULATIONSEIR_H_ */
