#include "RiemannianGraph.hpp"
#include "utils/linalg.hpp"

RiemannianGraph::RiemannianGraph(
    const std::vector<std::array<double, 3>> &vertices,
    const Octree<oNode> &octree, int k) {

  _max_edge = -std::numeric_limits<double>::infinity();
  _adj_list = std::vector<std::set<int>>(vertices.size());

  for (int i = 0; i < vertices.size(); i++) {

    auto v = vertices[i];
    std::vector<int> ni = octree.kNearestNeighbors(v, k);

    int last_idx = ni[ni.size() - 1];
    _max_edge = std::max(_max_edge, distance(v, vertices[last_idx]));

    for (int n : ni) {
      _adj_list[i].insert(n);
      _adj_list[n].insert(i);
    };
  };
};
