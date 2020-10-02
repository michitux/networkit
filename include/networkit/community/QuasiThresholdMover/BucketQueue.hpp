#ifndef BUCKETQUEUE_H
#define BUCKETQUEUE_H

#include <networkit/community/QuasiThresholdMover/DynamicForest.hpp>
#include <networkit/graph/Graph.hpp>

namespace NetworKit {

namespace QuasiThresholdMoving {
class BucketQueue {
public:
    BucketQueue(count n = 0);
    void fill(const std::vector<node> &elements, const DynamicForest &dynamicForest);
    void insertParent(node p);
    node next();
    bool empty() const;
    std::string printQueue() const;

private:
    count nextNode;
    count currentBucket;
    count n;
    std::vector<node> nodes;
    // points to first element in the bucket
    std::vector<count> border;
};
} // namespace QuasiThresholdMoving

} // namespace NetworKit

#endif // BUCKETQUEUE_H
