/*
 * GlobalCurveball.h
 *
 *  Created on: 26.05.2018
 *      Author: Manuel Penschuck <networkit@manuel.jetzt>
 */
#ifndef RANDOMIZATION_GLOBAL_CURVEBALL_H_
#define RANDOMIZATION_GLOBAL_CURVEBALL_H_

#include <memory>
#include <utility>

#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {

// pImpl
namespace CurveballDetails { struct GlobalCurveballImpl; }


class GlobalCurveball : public Algorithm {
public:
    /**
     * Instantiate a GlobalCurveball object
     *
     * @param G                        Undirected and unweighted graph to be randomized
     * @param number_of_global_trades  Number of global trades to be executed (each edge
     *                                 is considered exactly twice per global traded)
     * @param degreePreservingShufflePreprocessing Execute DegreePreservingShuffle
     *                                 (see Algorithm for description) as a preprocessing
     *                                 step. This is more efficient than calling the algorithm
     *                                 explicitly.
     */
    explicit GlobalCurveball(const Graph &G,
                             unsigned number_of_global_trades = 20,
                             bool degreePreservingShufflePreprocessing = false);


    virtual ~GlobalCurveball();

    /**
     * Execute trades as configured in the constructor.
     * @warning This function has to be called exactly one before invoking getGraph()
     */
    virtual void run() override final;

    /**
     * Returns a new graph instance with the same degree sequence as the input
     * graph, but with randomized neighbourhoods.
     */
    Graph getGraph();

    virtual std::string toString() const override final;

    virtual bool isParallel() const override final {
        return false;
    }

private:
    std::unique_ptr<CurveballDetails::GlobalCurveballImpl> impl;
    unsigned numGlobalTrades;
    bool degreePreservingShuffle;

};

}; // ! namespace NetworKit

#endif // ! RANDOMIZATION_GLOBAL_CURVEBALL_H_
