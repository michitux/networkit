/*
 * DotWriter.cpp
 *
 *  Created on: Jun 5, 2013
 *      Author: forigem
 */

#include <networkit/io/DotGraphWriter.hpp>

namespace NetworKit {

void DotGraphWriter::write(Graph& graph, std::string path) const {
	std::ofstream file{path};
	
	file << "graph {\n";
	graph.forEdges([&](node u, node v){
		file << u << " -- " << v << ";\n";
	});
	file << "}\n";
}

} /* namespace NetworKit */
