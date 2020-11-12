#ifndef QUASITHRESHOLDEDITINGLINEAR_H
#define QUASITHRESHOLDEDITINGLINEAR_H

#include <networkit/graph/Graph.hpp>

namespace NetworKit {
	namespace QuasiThresholdMoving {
		class QuasiThresholdEditingLinear {
		public:
			QuasiThresholdEditingLinear(const Graph& G);
			
			void run();
			
			std::vector<node> getParents() const;
			Graph getDefiningForest() const;
			Graph getQuasiThresholdGraph() const;
			
		private:
			const Graph& G;
			std::vector<node> parent;
			bool hasRun;
		};
		
	} // namespace NetworKit
	} //namespace QuasiThresholdMoving


#endif // QUASITHRESHOLDEDITINGLINEAR_H
