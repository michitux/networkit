#include "CKBDynamic.h"
#include "ckbdynamic/CKBDynamicImpl.h"


namespace NetworKit {
	CKBDynamic::CKBDynamic(count n, count minCommunitySize, count maxCommunitySize, double communitySizeExponent, double minSplitRatio, count minCommunityMembership, count maxCommunityMembership, double communityMembershipExponent, double communityEventProbability, double nodeEventProbability, double perturbationProbability, double intraCommunityEdgeProbability, double intraCommunityEdgeExponent, double epsilon, count numTimesteps) :
		parameters({n, minCommunitySize, maxCommunitySize, communitySizeExponent, minSplitRatio, minCommunityMembership, maxCommunityMembership, communityMembershipExponent, communityEventProbability, nodeEventProbability, perturbationProbability, intraCommunityEdgeProbability, intraCommunityEdgeExponent, epsilon, numTimesteps}) {}

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
