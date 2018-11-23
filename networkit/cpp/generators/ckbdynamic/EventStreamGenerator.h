#include "../../Globals.h"
#include "../../dynamics/GraphEvent.h"
#include "../../dynamics/CommunityEvent.h"
#include "../../base/Algorithm.h"

namespace NetworKit {
	namespace CKBDynamicImpl {
		class EventStreamGenerator : public Algorithm {
		public:
			EventStreamGenerator(count numTimesteps);

			void addEdge(index timestep, node u, node v);
			void removeEdge(index timestep, node u, node v);
			void addNode(index timestep, node u);
			void removeNode(index timestep, node u);
			void nodeJoinsCommunity(index timestep, node u, index community);
			void nodeLeavesCommunity(index timestep, node u, index community);

			virtual void run() override;

			std::vector<GraphEvent> getGraphEvents();
			std::vector<CommunityEvent> getCommunityEvents();
		private:
			struct InitialEvent {
				std::vector<std::pair<node, node>> edgesAdded;
				std::vector<std::pair<node, node>> edgesRemoved;
				std::vector<node> nodesAdded;
				std::vector<node> nodesRemoved;
				std::vector<std::pair<node, index>> nodeJoinsCommunity;
				std::vector<std::pair<node, index>> nodeLeavesCommunity;
			};

			std::vector<InitialEvent> initialEvents;

			// The time step when a node was born
			std::vector<index> birthTime;

			std::vector<GraphEvent> graphEvents;
			std::vector<CommunityEvent> communityEvents;
		};
	}
}
