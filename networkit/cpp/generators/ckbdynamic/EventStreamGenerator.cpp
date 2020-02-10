#include "EventStreamGenerator.h"
#include "NodePairHash.h"
#include <unordered_map>
#include "Community.h"


namespace NetworKit {
	namespace CKBDynamicImpl {
		EventStreamGenerator::EventStreamGenerator(count numTimesteps) : initialEvents(numTimesteps) {}

		void EventStreamGenerator::addEdge(index timestep, node u, node v) {
			initialEvents[timestep].edgesAdded.emplace_back(u, v);
		}

		void EventStreamGenerator::removeEdge(index timestep, node u, node v) {
			initialEvents[timestep].edgesRemoved.emplace_back(u, v);
		}

		void EventStreamGenerator::addNode(index timestep, node u) {
			assert(u == birthTime.size());
			birthTime.push_back(timestep);
			initialEvents[timestep].nodesAdded.emplace_back(u);
		}

		void EventStreamGenerator::removeNode(index timestep, node u) {
			initialEvents[timestep].nodesRemoved.emplace_back(u);
		}

		void EventStreamGenerator::nodeJoinsCommunity(index timestep, node u, index community) {
			initialEvents[timestep].nodeJoinsCommunity.emplace_back(u, community);
		}

		void EventStreamGenerator::nodeLeavesCommunity(index timestep, node u, index community) {
			initialEvents[timestep].nodeLeavesCommunity.emplace_back(u, community);
		}

		std::vector<GraphEvent> EventStreamGenerator::getGraphEvents() {
			assureFinished();

			return std::move(graphEvents);
		}

		std::vector<CommunityEvent> EventStreamGenerator::getCommunityEvents() {
			assureFinished();

			return std::move(communityEvents);
		}

		void EventStreamGenerator::run() {
			std::vector<std::unordered_map<node, count>> neighbors(birthTime.size());
			std::vector<index> deathTime(birthTime.size(), none);

			std::unordered_set<std::pair<node, node>, NodePairHash> freshlyAddedEdges;
			std::unordered_map<std::pair<index, node>, count, NodePairHash> freshlyJoinedCommunities;

			std::vector<std::pair<node, node>> edgesAdded, edgesRemoved;
			for (index timestep = 0; timestep < initialEvents.size(); ++timestep) {
				for (std::pair<node, node> e : initialEvents[timestep].edgesAdded) {
					index maxBirth = std::max(birthTime[e.first], birthTime[e.second]);
					if (maxBirth > timestep) {
						// Delay edge creation to when both nodes exist
						initialEvents[maxBirth].edgesAdded.push_back(e);
					} else {
						auto it = neighbors[e.first].find(e.second);
						if (it == neighbors[e.first].end()) {
							neighbors[e.first].emplace(e.second, 1);
							neighbors[e.second].emplace(e.first, 1);
							freshlyAddedEdges.insert(e);
						} else {
							++(it->second);
							it = neighbors[e.second].find(e.first);
							assert(it != neighbors[e.second].end());
							++(it->second);
						}
					}
				}

				for (std::pair<node, node> e : initialEvents[timestep].edgesRemoved) {
					index minDeath = std::min(deathTime[e.first], deathTime[e.second]);
					if (minDeath < timestep) {
						assert(neighbors[e.first].find(e.second) == neighbors[e.first].end());
						assert(neighbors[e.second].find(e.first) == neighbors[e.second].end());
					} else {
						auto it = neighbors[e.first].find(e.second);
						assert(it != neighbors[e.first].end());
						--(it->second);

						if (it->second == 0) {
							neighbors[e.first].erase(it);
							auto eit = freshlyAddedEdges.find(e);
							if (eit != freshlyAddedEdges.end()) {
								freshlyAddedEdges.erase(eit);
							} else {
								graphEvents.emplace_back(GraphEvent::EDGE_REMOVAL, e.first, e.second);
							}
						}

						it = neighbors[e.second].find(e.first);
						assert(it != neighbors[e.second].end());
						--(it->second);

						if (it->second == 0) {
							assert(neighbors[e.first].find(e.second) == neighbors[e.first].end());
							neighbors[e.second].erase(it);
						} else {
							assert(it->second == neighbors[e.first][e.second]);
						}
					}
				}

				for (node u : initialEvents[timestep].nodesRemoved) {
					for (auto it : neighbors[u]) {
						node v = it.first;
						neighbors[v].erase(u);
						std::pair<node, node> e = Community::canonicalEdge(u, v);
						auto eit = freshlyAddedEdges.find(e);
						if (eit != freshlyAddedEdges.end()) {
							freshlyAddedEdges.erase(eit);
						} else {
							graphEvents.emplace_back(GraphEvent::EDGE_REMOVAL, e.first, e.second);
						}
					}

					neighbors[u].clear();
					deathTime[u] = timestep;

					graphEvents.emplace_back(GraphEvent::NODE_REMOVAL, u);
				}

				for (node u : initialEvents[timestep].nodesAdded) {
					graphEvents.emplace_back(GraphEvent::NODE_ADDITION, u);
				}

				for (auto it : freshlyAddedEdges) {
					graphEvents.emplace_back(GraphEvent::EDGE_ADDITION, it.first, it.second);
				}

				freshlyAddedEdges.clear();

				for (std::pair<node, index> nodeCom : initialEvents[timestep].nodeJoinsCommunity) {
					// Nodes might be added several times to the same community.
					// A possible reason is that during a split event, first the community is duplicated and thus a node joins it.
					// Then, it is removed from the community to remove overlap.
					// Later, it is added again because more nodes shall be added to the community.
					++freshlyJoinedCommunities[nodeCom];
				}

				for (std::pair<node, index> nodeCom : initialEvents[timestep].nodeLeavesCommunity) {
					auto eit = freshlyJoinedCommunities.find(nodeCom);
					if (eit != freshlyJoinedCommunities.end()) {
						if (eit->second > 1) {
							eit->second--;
						} else {
							freshlyJoinedCommunities.erase(eit);
						}
					} else {
						communityEvents.emplace_back(CommunityEvent::NODE_LEAVES_COMMUNITY, nodeCom.first, nodeCom.second);
					}
				}

				for (auto it : freshlyJoinedCommunities) {
					assert(it.second == 1);
					communityEvents.emplace_back(CommunityEvent::NODE_JOINS_COMMUNITY, it.first.first, it.first.second);
				}

				freshlyJoinedCommunities.clear();

				graphEvents.emplace_back(GraphEvent::TIME_STEP);
				communityEvents.emplace_back(CommunityEvent::TIME_STEP);
			}

			hasRun = true;
		}

	}
}
