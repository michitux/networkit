#include "CommunityMergeEvent.h"
#include "CKBDynamicImpl.h"
#include "../../auxiliary/Random.h"
#include <tlx/unused.hpp>

namespace NetworKit {
	namespace CKBDynamicImpl {
		CommunityMergeEvent::CommunityMergeEvent(CommunityPtr communityA, CommunityPtr communityB, count numSteps, CKBDynamicImpl& generator) : CommunityChangeEvent(generator, numSteps), targetEdgeProbability(.0), communityA(communityA), communityB(communityB) {
			count overlapSize = 0;

			for (node u : communityA->getNodes()) {
				if (communityB->hasNode(u)) {
					++overlapSize;
				} else {
					nodesToAddToB.push_back(u);
				}
			}

			std::shuffle(nodesToAddToB.begin(), nodesToAddToB.end(), Aux::Random::getURNG());

			for (node u : communityB->getNodes()) {
				if (!communityA->hasNode(u)) {
					nodesToAddToA.push_back(u);
				}
			}

			std::shuffle(nodesToAddToA.begin(), nodesToAddToA.end(), Aux::Random::getURNG());

			std::tie(std::ignore, targetEdgeProbability) = generator.communitySizeSampler->mergeCommunities(communityA->getNumberOfNodes(), communityA->getEdgeProbability(), communityB->getNumberOfEdges(), communityB->getEdgeProbability(), communityA->getNumberOfNodes() + communityB->getNumberOfNodes() - overlapSize);
			targetEdgeProbabilityPerCommunity = 1 - std::sqrt(1 - targetEdgeProbability);

			communityA->setAvailable(false);
			communityB->setAvailable(false);
		}

		void CommunityMergeEvent::nextStep() {
			auto addNodes = [&](CommunityPtr com, std::vector<node>& nodes) {
						count numNodesToAdd = nodes.size() / (numSteps - currentStep);
						for (index i = 0; i < numNodesToAdd; ++i) {
							node u = nodes.back();
							nodes.pop_back();

							// skip nodes that have been removed in the meantime
							if (generator.hasNode(u)) {
								com->addNode(u);
							}
						}
					};

			addNodes(communityA, nodesToAddToA);
			addNodes(communityB, nodesToAddToB);

			adaptProbability(communityA, targetEdgeProbabilityPerCommunity);
			adaptProbability(communityB, targetEdgeProbabilityPerCommunity);

			++currentStep;
			if (currentStep == numSteps) {
				active = false;
				const count oldNodes = communityA->getNumberOfNodes();
				tlx::unused(oldNodes);
				communityA->combineWith(*communityB);
				generator.removeCommunity(communityB);
				assert(communityA->getNumberOfNodes() == oldNodes);
				// This shouldn't change much but otherwise the community will loose half of
				// the edges in the next perturbation.
				communityA->changeEdgeProbability(targetEdgeProbability);
				communityA->setAvailable(true);
			}
		}
	}
}
