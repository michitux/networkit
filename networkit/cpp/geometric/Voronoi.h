#ifndef VORONOI_H_
#define VORONOI_H_

#ifndef NOOPENCV
#include <opencv2/opencv.hpp>
#endif
#include "Point2D.h"
#include "../base/Algorithm.h"
#include <libcola/cola.h>

namespace NetworKit {

  class Voronoi : public Algorithm {
  public:
    Voronoi(const std::vector<std::pair<float, float>>& points, float xmin, float xmax, float ymin, float ymax);
    virtual void run() override;

    std::vector<std::vector<std::pair<float, float>>> getFacets();

  private:
#ifndef NOOPENCV
    cv::Subdiv2D subdivision;
#endif
  };

} // namespace NetworKit

#endif
