#include "CKBDynamic.h"
#include "ckbdynamic/CKBDynamicImpl.h"


namespace NetworKit {
	CKBDynamic::CKBDynamic(count n, count minCommunitySize, count maxCommunitySize, double communitySizeExponent, count minCommunityMembership, count maxCommunityMembership, double communityMembershipExponent, double communityEventProbability, double nodeEventProbability, double perturbationProbability, double intraCommunityEdgeProbability, double intraCommunityEdgeExponent, double epsilon, double edgeSharpness, count numTimesteps) :
		parameters({n, minCommunitySize, maxCommunitySize, communitySizeExponent, minCommunityMembership, maxCommunityMembership, communityMembershipExponent, communityEventProbability, nodeEventProbability, perturbationProbability, intraCommunityEdgeProbability, intraCommunityEdgeExponent, epsilon, edgeSharpness, numTimesteps, nullptr, nullptr}) {}

	CKBDynamic::CKBDynamic(count n, const Graph& G, const Cover& C, double communityEventProbability, double nodeEventProbability, double perturbationProbability, double edgeSharpness, count numTimesteps) : parameters({n, 0, 0, .0, 0, 0, .0, communityEventProbability, nodeEventProbability, perturbationProbability, .0, .0, .0, edgeSharpness, numTimesteps, &G, &C}) {}

	void CKBDynamic::run() {
		hasRun = false;
		CKBDynamicImpl::CKBDynamicImpl impl(parameters);
		impl.run();
		graphEvents = impl.getGraphEvents();
		communityEvents = impl.getCommunityEvents();
		hasRun = true;
	}


	const std::vector<GraphEvent>& CKBDynamic::getGraphEvents() const {
		assureFinished();
		return graphEvents;
	}

	const std::vector<CommunityEvent>& CKBDynamic::getCommunityEvents() const {
		assureFinished();
		return communityEvents;
	}
}
