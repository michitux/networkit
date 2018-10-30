#include "CommunityMergeEvent.h"
#include "CKBDynamicImpl.h"
#include "../../auxiliary/Random.h"
#include <tlx/unused.hpp>

namespace NetworKit {
	namespace CKBDynamicImpl {
		CommunityMergeEvent::CommunityMergeEvent(CommunityPtr communityA, CommunityPtr communityB, count targetSize, double targetEdgeProbability, count numSteps, CKBDynamicImpl& generator) : CommunityChangeEvent(generator, numSteps), communities({communityA, communityB}), targetSize(targetSize), targetEdgeProbability(targetEdgeProbability) {
			targetEdgeProbabilityPerCommunity = 1 - std::sqrt(1 - targetEdgeProbability);

			for (count c = 0; c < 2; ++c) {
				// Get a list of nodes to be added to community c
				for (node u : communities[1-c]->getNodes()) {
					if (!communities[c]->hasNode(u)) {
						nodesToAddTo[c].push_back(u);
					}
				}

				std::shuffle(nodesToAddTo[c].begin(), nodesToAddTo[c].end(), Aux::Random::getURNG());

				communities[c]->setAvailable(false);
			}

		}

		void CommunityMergeEvent::nextStep() {

			const count estimatedOverlap = communities[0]->getNumberOfNodes() - nodesToAddTo[1].size();
			const count estimatedOverlap1 = communities[1]->getNumberOfNodes() - nodesToAddTo[0].size();

			// FIXME: fix this for node deletions
			assert(estimatedOverlap == estimatedOverlap1);
			const count estimatedMergedSize = communities[0]->getNumberOfNodes() + nodesToAddTo[0].size();

			const count totalNodesToAdd = estimatedOverlap > targetSize ? 0 : (targetSize - estimatedOverlap) / (numSteps - currentStep);
			const count totalNodesToRemove = estimatedMergedSize < targetSize ? 0 : (estimatedMergedSize - targetSize) / (numSteps - currentStep);


			auto nonOverlappingNodes = [&]() -> count {
							   return nodesToAddTo[0].size() + nodesToAddTo[1].size();
						   };

			count nodesRemoved = 0;

			if (nonOverlappingNodes() > 0 && totalNodesToRemove > 0) {
				double fractionOfNodesToRemove = totalNodesToRemove * 1.0 / nonOverlappingNodes();
				for (count c = 0; c < 2; ++c) {
					count nodesToRemove = nodesToAddTo[c-1].size();

					if (fractionOfNodesToRemove < 1) {
						nodesToRemove = fractionOfNodesToRemove * nodesToAddTo[1-c].size();
						if (Aux::Random::real() > fractionOfNodesToRemove - static_cast<count>(fractionOfNodesToRemove)) {
							nodesToRemove += 1;
						}

						nodesToRemove = std::min(totalNodesToRemove - nodesRemoved, nodesToRemove);
					}

					for (count i = 0; i < nodesToRemove; ++i) {
						node u = nodesToAddTo[1-c].back();
						nodesToAddTo[1-c].pop_back();

						if (communities[c]->hasNode(u)) {
							communities[c]->removeNode(u);
						}

						++nodesRemoved;
					}
				}
			}

			if (nodesRemoved < totalNodesToRemove && nonOverlappingNodes() == 0) {
				for (; nodesRemoved < totalNodesToRemove; ++nodesRemoved) {
					assert(communities[0]->getNumberOfNodes() > 0);
					assert(communities[0]->getNumberOfNodes() == communities[1]->getNumberOfNodes());
					const node u = communities[0]->removeRandomNode();
					communities[1]->removeNode(u);
				}
			}


			for (count c = 0; c < 2; ++c) {
				// first adapt the probability so new
				// nodes get directly the right amount
				// of neighbors
				adaptProbability(communities[c], targetEdgeProbabilityPerCommunity);
			}

			count nodesAdded = 0;

			if (nonOverlappingNodes() > 0 && totalNodesToAdd > 0) {
				double fractionOfNodesToAdd = totalNodesToAdd * 1.0 / nonOverlappingNodes();

				for (count c = 0; c < 2; ++c) {
					count numNodesToAdd = nodesToAddTo[c].size();

					if (fractionOfNodesToAdd < 1) {
						numNodesToAdd = fractionOfNodesToAdd * nodesToAddTo[c].size();
						if (Aux::Random::real() > fractionOfNodesToAdd - static_cast<count>(fractionOfNodesToAdd)) {
							numNodesToAdd += 1;
						}

						numNodesToAdd = std::min(totalNodesToAdd - nodesAdded, numNodesToAdd);
					}

					for (count i = 0; i < numNodesToAdd; ++i) {
						node u = nodesToAddTo[c].back();
						nodesToAddTo[c].pop_back();

						// skip nodes that have been removed in the meantime
						if (generator.hasNode(u)) {
							communities[c]->addNode(u);
						}

						++nodesAdded;
					}
				}
			}

			if (nodesAdded < totalNodesToAdd && nonOverlappingNodes() == 0) {
				count numNodesToAdd = totalNodesToAdd - nodesAdded;
				std::vector<node> nodes = generator.communityNodeSampler.birthCommunityNodes(numNodesToAdd, communities[0]->getNodes());
				for (node u : nodes) {
					communities[0]->addNode(u);
					communities[1]->addNode(u);
				}
			}

			++currentStep;
			if (currentStep == numSteps) {
				active = false;
				const count oldNodes = communities[0]->getNumberOfNodes();
				tlx::unused(oldNodes);
				communities[0]->combineWith(*communities[1]);
				generator.removeCommunity(communities[1]);
				assert(communities[0]->getNumberOfNodes() == oldNodes);
				// This shouldn't change much but otherwise the community will loose half of
				// the edges in the next perturbation.
				communities[0]->changeEdgeProbability(targetEdgeProbability);
				communities[0]->setAvailable(true);
			}
		}
	}
}
