/*
 *
 */

#include <networkit/community/QuasiThresholdMover/BucketQueue.hpp>

namespace NetworKit {

namespace QuasiThresholdMoving {
BucketQueue::BucketQueue(count n)
    : n(n), nodes(n), border(n), nextNode(none), currentBucket(none) {}

void BucketQueue::fill(const std::vector<node> &elements, const DynamicForest &dynamicForest) {
    count maxDepth = std::min(n - 1, 2 * elements.size());
    auto borderEnd = border.begin() + std::min(maxDepth + 1, border.size());
    std::fill(border.begin(), borderEnd, 0);
    for (node u : elements) {
        if (dynamicForest.depth(u) > maxDepth)
            continue;
        border[dynamicForest.depth(u)] += 1;
    }
    std::partial_sum(border.begin(), borderEnd, border.begin());
    currentBucket = maxDepth;
    nextNode = none;
    std::vector<node>::const_reverse_iterator rit = elements.rbegin();
    for (; rit != elements.rend(); ++rit) {
        node u = *rit;
        if (dynamicForest.depth(u) > maxDepth) {
            continue;
        }
        border[dynamicForest.depth(u)] -= 1;
        nodes[border[dynamicForest.depth(u)]] = u;
        nextNode += 1;
    }
}

node BucketQueue::next() {
    if (nextNode == none) {
        return none;
    }
    node result = nodes[nextNode];
    while (nextNode < border[currentBucket]) {
        currentBucket -= 1;
    }
    nextNode -= 1;
    return result;
}

void BucketQueue::insertParent(node p) {
    nextNode += 1;
    // first element of currentBucket
    count bucketBorder = border[currentBucket];
    node firstOfBucket = nodes[bucketBorder];
    nodes[nextNode] = firstOfBucket;
    nodes[bucketBorder] = p;
    border[currentBucket] += 1;
}

bool BucketQueue::empty() const {
    return nextNode == none;
}

std::string BucketQueue::printQueue() const {
    if (empty()) {
        return "BucketQueue:";
    }
    std::stringstream ss;
    ss << "BucketQueue:";
    int j = 0;
    for (count i = 0; i <= nextNode; i++) {
        while (border[j] == i) {
            ss << "| ";
            j++;
        }
        ss << nodes[i] << " ";
    }
    return ss.str();
}
} // namespace QuasiThresholdMoving

} // namespace NetworKit
