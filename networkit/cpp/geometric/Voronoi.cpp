#include "Voronoi.h"


namespace NetworKit {

#ifndef NOOPENCV
  Voronoi::Voronoi(const std::vector<std::pair<float, float>>& points, float xmin, float xmax, float ymin, float ymax) : subdivision(cv::Rect2f(xmin, ymin, xmax-xmin, ymax-ymin)) {
    for (auto p : points) {
      cv::Point2f p2f(p.first, p.second);
      subdivision.insert(p2f);
    }
  }
#else
  Voronoi::Voronoi(const std::vector<std::pair<float, float>>& points, float xmin, float xmax, float ymin, float ymax) {
    throw std::runtime_error("NetworKit compiled without OpenCV, no Voronoi calculation available.");
  }
#endif

  void Voronoi::run() {
    // Nothing to do?!
    hasRun = true;
  }


  std::vector<std::vector<std::pair<float, float>>> Voronoi::getFacets() {
#ifndef NOOPENCV
    std::vector<std::vector<cv::Point2f>> facetList;
    std::vector<cv::Point2f> centerList;
    subdivision.getVoronoiFacetList(std::vector<int>(), facetList, centerList);

    std::vector<std::vector<std::pair<float, float>>> result;
    result.resize(facetList.size());

    for (size_t i = 0; i < facetList.size(); ++i) {
      result[i].reserve(facetList[i].size());
      for (cv::Point2f p : facetList[i]) {
	result[i].push_back(std::make_pair(p.x, p.y));
      }
    }

    return result;
#else
    // unreachable as the constructor throws!
    die();
#endif
  }
}
