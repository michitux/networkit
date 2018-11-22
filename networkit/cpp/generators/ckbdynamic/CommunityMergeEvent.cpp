#include "CommunityMergeEvent.h"
#include "CKBDynamicImpl.h"
#include "../../auxiliary/Random.h"
#include <tlx/unused.hpp>

namespace NetworKit {
	namespace CKBDynamicImpl {
		CommunityMergeEvent::CommunityMergeEvent(CommunityPtr communityA, CommunityPtr communityB, count targetSize, double targetEdgeProbability, count numSteps, CKBDynamicImpl& generator) : CommunityChangeEvent(generator, numSteps), communities({communityA, communityB}), targetSize(targetSize), targetEdgeProbability(targetEdgeProbability), communitiesMerged(false) {
			targetEdgeProbabilityPerCommunity = 1 - std::sqrt(1 - targetEdgeProbability);

			for (count c = 0; c < 2; ++c) {
				communities[c]->setCurrentEvent(this);

				// Get a list of nodes to be added to community c
				for (node u : communities[1-c]->getNodes()) {
					if (!communities[c]->hasNode(u)) {
						nodesToAddTo[c].insert(u);
					}
				}
			}

		}

		void CommunityMergeEvent::nextStep() {

			const count estimatedOverlap = communities[0]->getNumberOfNodes() - nodesToAddTo[1].size();

			assert(communitiesMerged || estimatedOverlap == communities[1]->getNumberOfNodes() - nodesToAddTo[0].size());
			const count estimatedMergedSize = communities[0]->getNumberOfNodes() + nodesToAddTo[0].size();
			assert(communitiesMerged || estimatedMergedSize == communities[1]->getNumberOfNodes() + nodesToAddTo[1].size());

			const count totalNodesToAdd = estimatedOverlap > targetSize ? 0 : (targetSize - estimatedOverlap) / (numSteps - currentStep);
			const count totalNodesToRemove = estimatedMergedSize < targetSize ? 0 : (estimatedMergedSize - targetSize) / (numSteps - currentStep);


			auto nonOverlappingNodes = [&]() -> count {
							   return nodesToAddTo[0].size() + nodesToAddTo[1].size();
						   };

			count nodesRemoved = 0;

			// Reset the desired number of nodes to the actual number of nodes.
			// If any change shall happen, it will be requested below.
			communities[0]->setDesiredNumberOfNodes(communities[0]->getNumberOfNodes());
			if (!communitiesMerged) {
				communities[1]->setDesiredNumberOfNodes(communities[1]->getNumberOfNodes());
			}

			if (nonOverlappingNodes() > 0 && totalNodesToRemove > 0) {
				std::array<count, 2> nodesToRemove {nodesToAddTo[0].size(), nodesToAddTo[1].size()};
				if (nodesToRemove[0] + nodesToRemove[1] > totalNodesToRemove) {
					nodesToRemove[0] = nodesToAddTo[0].size() * totalNodesToRemove / nonOverlappingNodes();
					nodesToRemove[1] = totalNodesToRemove - nodesToRemove[0];

					if (nodesToRemove[1] > nodesToAddTo[1].size() || (nodesToAddTo[0].size() > nodesToRemove[0] && nodesToRemove[1] > 0 && Aux::Random::real() < (nodesToAddTo[0].size() *  totalNodesToRemove * 1.0 / nonOverlappingNodes()) - nodesToRemove[0])) {
						++nodesToRemove[0];
						--nodesToRemove[1];
					}

					assert(nodesToRemove[0] + nodesToRemove[1] == totalNodesToRemove);
				}

				assert(nodesToRemove[0] <= nodesToAddTo[0].size());
				assert(nodesToRemove[1] <= nodesToAddTo[1].size());

				for (count c = 0; c < 2; ++c) {
					for (count i = 0; i < nodesToRemove[c]; ++i) {
						const node u = nodesToAddTo[c].at(Aux::Random::index(nodesToAddTo[c].size()));
						nodesToAddTo[c].erase(u);

						communities[1-c]->removeNode(u);
						++nodesRemoved;
					}

					communities[1-c]->setDesiredNumberOfNodes(communities[1-c]->getNumberOfNodes());
				}
			}

			if (nonOverlappingNodes() == 0) {
				mergeCommunities();
			}

			if (nodesRemoved < totalNodesToRemove) {
				assert(nonOverlappingNodes() == 0);
				for (; nodesRemoved < totalNodesToRemove; ++nodesRemoved) {
					assert(communities[0]->getNumberOfNodes() > 0);
					// communities[1] doesn't exist anymore
					communities[0]->removeRandomNode();
				}

				communities[0]->setDesiredNumberOfNodes(communities[0]->getNumberOfNodes());
			}


			if (communitiesMerged) {
				adaptProbability(communities[0], targetEdgeProbability);
			} else {
				for (count c = 0; c < 2; ++c) {
					// first adapt the probability so new
					// nodes get directly the right amount
					// of neighbors
					adaptProbability(communities[c], targetEdgeProbabilityPerCommunity);
				}
			}

			count nodesAdded = 0;

			// At minimum, we need to add enough nodes such that both communities reach the minimum size.
			const count minSize = generator.communitySizeSampler->getMinSize();
			std::array<count, 2> minNodesToAdd {0,0};
			for (count c = 0; c < 2; ++c) {
				if (communities[c]->getNumberOfNodes() < minSize) {
					minNodesToAdd[c] = minSize - communities[c]->getNumberOfNodes();
				}
				if (communitiesMerged) break;
			}

			assert(communitiesMerged || communities[0]->getNumberOfNodes() + nodesToAddTo[0].size() == communities[1]->getNumberOfNodes() + nodesToAddTo[1].size());

			if (nonOverlappingNodes() > 0 && (totalNodesToAdd > 0 || minNodesToAdd[0] + minNodesToAdd[1] > 0)) {
				std::array<count, 2> numNodesToAdd {nodesToAddTo[0].size(), nodesToAddTo[1].size()};

				// If we have more non-overlapping nodes than nodes to be added, select some.
				// However, select only adding all non-overlapping nodes achieves the minimum size.
				if (numNodesToAdd[0] + numNodesToAdd[1] > totalNodesToAdd && communities[0]->getNumberOfNodes() + nodesToAddTo[0].size() > minSize) {
					numNodesToAdd[0] = numNodesToAdd[0] * totalNodesToAdd / nonOverlappingNodes();
					numNodesToAdd[1] = totalNodesToAdd - numNodesToAdd[0];

					// Always round the other way if the previous allocation is infeasible.
					// Apply probabilitic rounding if there is enough possibility for change.
					if (numNodesToAdd[1] > nodesToAddTo[1].size() || (nodesToAddTo[0].size() > numNodesToAdd[0] && numNodesToAdd[1] > 0 && Aux::Random::real() < nodesToAddTo[0].size() *  totalNodesToAdd *  1.0 / nonOverlappingNodes() - numNodesToAdd[0])) {
						++numNodesToAdd[0];
						--numNodesToAdd[1];
					}

					// Ensure that if possible, we are adding at least the minimum amount of nodes required to reach the minimum community size.
					for (count c = 0; c < 2; ++c) {
						if (numNodesToAdd[c] < minNodesToAdd[c]) {
							numNodesToAdd[c] = minNodesToAdd[c];
							numNodesToAdd[1-c] = std::max(totalNodesToAdd > numNodesToAdd[c] ? totalNodesToAdd - numNodesToAdd[c] : 0, minNodesToAdd[1 - c]);
						}
					}
					assert(numNodesToAdd[0] + numNodesToAdd[1] >= totalNodesToAdd);
					assert(numNodesToAdd[0] <= nodesToAddTo[0].size());
					assert(numNodesToAdd[1] <= nodesToAddTo[1].size());
				}

				for (count c = 0; c < 2; ++c) {
					for (count i = 0; i < numNodesToAdd[c]; ++i) {
						const node u = nodesToAddTo[c].at(Aux::Random::index(nodesToAddTo[c].size()));

						assert(communities[1-c]->hasNode(u));
						communities[c]->addNode(u);
						nodesToAddTo[c].erase(u);
						++nodesAdded;
					}

					communities[c]->setDesiredNumberOfNodes(communities[c]->getNumberOfNodes());
				}
			}

			if (nodesAdded < totalNodesToAdd || communities[0]->getNumberOfNodes() < minSize) {
				// Assert the overlap is empty now
				if (nonOverlappingNodes() > 0) throw std::logic_error("There are nodes not in the overlap but we did not add them even though we should");

				mergeCommunities();

				count numNodesToAdd = totalNodesToAdd - nodesAdded;
				if (communities[0]->getNumberOfNodes() + numNodesToAdd < minSize) {
					numNodesToAdd = minSize - communities[0]->getNumberOfNodes();
				}

				communities[0]->setDesiredNumberOfNodes(communities[0]->getNumberOfNodes() + numNodesToAdd);
			}

			++currentStep;
			if (currentStep == numSteps) {
				mergeCommunities();

				communities[0]->changeEdgeProbability(targetEdgeProbability);
				communities[0]->setCurrentEvent(nullptr);
				active = false;
			}
		}

		void CommunityMergeEvent::mergeCommunities() {
			if (communitiesMerged) return;

			communities[1]->setCurrentEvent(nullptr);

			const count oldNodes = communities[0]->getNumberOfNodes();
			tlx::unused(oldNodes);
			communities[0]->combineWith(*communities[1]);
			generator.removeCommunity(communities[1]);
			assert(communities[0]->getNumberOfNodes() == oldNodes);

			// Set the edge probability to the probability that should match the combined number of edges
			double probNonEdge = 1.0 - communities[0]->getEdgeProbability();
			communities[0]->changeEdgeProbability(1.0 - probNonEdge*probNonEdge);

			communitiesMerged = true;
		}

		void CommunityMergeEvent::notifyNodeRemovedFromCommunity(node u, CommunityPtr com) {
			assert(com == communities[0] || com == communities[1]);

			for (count c = 0; c < 2; ++c) {
				nodesToAddTo[c].erase(u);
			}
		}

		void CommunityMergeEvent::notifyNodeAddedToCommunity(node u, CommunityPtr com) {
			if (!communitiesMerged) {
				if (com == communities[0]) {
					assert(nodesToAddTo[0].contains(u));
				} else {
					assert(nodesToAddTo[1].contains(u));
				}
			}
		}

		bool CommunityMergeEvent::canRemoveNode() const {
			return communitiesMerged;
		}
	}
}
