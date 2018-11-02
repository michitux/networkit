#ifndef CKB_DYNAMIC_IMPL_H
#define CKB_DYNAMIC_IMPL_H

#include "../../graph/Graph.h"
#include "../../structures/Cover.h"
#include "../../dynamics/GraphEvent.h"
#include "../../dynamics/CommunityEvent.h"
#include "../../base/Algorithm.h"
#include "../../auxiliary/SamplingSet.h"
#include "Community.h"
#include "BucketSampling.h"
#include "CommunitySizeDistribution.h"
#include "CommunityChangeEvent.h"
#include "../CKBDynamic.h"
#include "NodePairHash.h"

namespace NetworKit {
	namespace CKBDynamicImpl {
		class CKBDynamicImpl : public Algorithm {
		public:
			CKBDynamicImpl(const Graph& G, const Cover& model, count numTimesteps);
			CKBDynamicImpl(const CKBDynamic::param_type& parameters);

			virtual void run() override;

			std::vector<GraphEvent> getGraphEvents() const;
			std::vector<CommunityEvent> getCommunityEvents() const;
		private:
			std::vector<GraphEvent> graphEvents;
			std::vector<CommunityEvent> communityEvents;
		public:

			void addEdge(node u, node v);
			void removeEdge(node u, node v);
			void addNodeToCommunity(node u, CommunityPtr com);
			void removeNodeFromCommunity(node u, CommunityPtr com);
			void addCommunity(CommunityPtr com);
			void removeCommunity(CommunityPtr com);

			void generateNode();
			void eraseNode();

			index nextCommunityId();
			count sampleNumSteps() const;

			bool hasNode(node u) const { return nodesAlive.contains(u); };

			BucketSampling communityNodeSampler;
			std::unique_ptr<CommunitySizeDistribution> communitySizeSampler;
		private:
			void finishTimeStep();

			Aux::SamplingSet<CommunityPtr> availableCommunities;
			Aux::SamplingSet<CommunityPtr> communities;
			std::vector<std::unordered_set<CommunityPtr>> nodeCommunities;
			CommunityPtr globalCommunity;
			count maxCommunityId;

			Aux::SamplingSet<node> nodesAlive;
			std::unordered_map<std::pair<node, node>, count, NodePairHash> edgesAlive;

			// Events of the current time step
			std::unordered_map<std::pair<node, node>, GraphEvent, NodePairHash> currentEdgeEvents;
			std::vector<node> currentErasedNodes;
			std::unordered_map<std::pair<index, node>, CommunityEvent, NodePairHash> currentCommunityEvents;

			count n;
			double communityEventProbability;
			double nodeEventProbability;
			double perturbationProbability;
			double epsilon;
			count numTimesteps;
			count currentCommunityMemberships;
			std::vector<std::unique_ptr<CommunityChangeEvent>> currentEvents;
		};

	}
} // namespace NetworKit

#endif
