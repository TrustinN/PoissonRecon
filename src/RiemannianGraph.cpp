#include "RiemannianGraph.hpp"
#include "utils.hpp"

RiemannianGraph::RiemannianGraph(
    const std::vector<std::array<double, 3>> &vertices, const Octree &octree,
    int k) {
  max_edge = -std::numeric_limits<double>::infinity();

  for (int i = 0; i < vertices.size(); i++) {
    auto v = vertices[i];
    std::vector<int> ni = octree.kNearestNeighbors(v, k);
    std::vector<std::array<double, 3>> np;
    for (int n : ni) {
      np.push_back(vertices[n]);
    }

    max_edge = std::max(max_edge, distance(v, np[np.size() - 1]));
    for (int n : ni) {
      adj_list[i].insert(n);
      adj_list[n].insert(i);
    };
  };
};
