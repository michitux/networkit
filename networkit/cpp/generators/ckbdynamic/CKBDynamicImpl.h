#ifndef CKB_DYNAMIC_IMPL_H
#define CKB_DYNAMIC_IMPL_H

#include "../../graph/Graph.h"
#include "../../structures/Cover.h"
#include "../../dynamics/GraphEvent.h"
#include "../../dynamics/CommunityEvent.h"
#include "../../base/Algorithm.h"
#include "../../auxiliary/SamplingSet.h"
#include "Community.h"
#include "CommunitySizeDistribution.h"
#include "CommunityChangeEvent.h"
#include "../CKBDynamic.h"
#include "NodePairHash.h"
#include "EventStreamGenerator.h"
#include "CommunityMembershipDistribution.h"

namespace NetworKit {
	namespace CKBDynamicImpl {
		class CKBDynamicImpl : public Algorithm {
		public:
			CKBDynamicImpl(const Graph& G, const Cover& model, count numTimesteps);
			CKBDynamicImpl(const CKBDynamic::param_type& parameters);

			virtual void run() override;

			std::vector<GraphEvent> getGraphEvents();
			std::vector<CommunityEvent> getCommunityEvents();

			void addEdge(node u, node v, bool nodeJoined);
			void removeEdge(node u, node v, bool nodeLeft);
			void addNodeToCommunity(node u, CommunityPtr com);
			void removeNodeFromCommunity(node u, CommunityPtr com);
			void addCommunity(CommunityPtr com);
			void removeCommunity(CommunityPtr com);

			void generateNode();
			void eraseNode();

			index nextCommunityId();
			count sampleNumSteps() const;

			bool hasNode(node u) const { return nodesAlive.contains(u); };

			std::unique_ptr<CommunitySizeDistribution> communitySizeSampler;
		private:
			std::unique_ptr<CommunityMembershipDistribution> membershipDistribution;
			void assignNodesToCommunities();

			Aux::SamplingSet<CommunityPtr> availableCommunities;
			Aux::SamplingSet<CommunityPtr> communities;
			std::vector<std::unordered_set<CommunityPtr>> nodeCommunities;
			CommunityPtr globalCommunity;
			count maxCommunityId;

			Aux::SamplingSet<node> nodesAlive;
			std::vector<count> desiredMemberships;
			count sumOfDesiredMemberships;

			count currentTimeStep;
			EventStreamGenerator eventStream;

			count n;
			double communityEventProbability;
			double nodeEventProbability;
			double perturbationProbability;
			double epsilon;
			double edgeSharpness;
			count numTimesteps;
			count currentCommunityMemberships;
			std::vector<std::unique_ptr<CommunityChangeEvent>> currentEvents;
		};

	}
} // namespace NetworKit

#endif
