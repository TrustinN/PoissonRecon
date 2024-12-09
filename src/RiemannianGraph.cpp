#include "RiemannianGraph.hpp"
#include <omp.h>

RiemannianGraph::RiemannianGraph(
    const std::vector<std::array<double, 3>> &vertices, const Octree &octree,
    int k) {

  int num_vertices = vertices.size();
  _adj_list = std::vector<std::set<int>>(num_vertices);
  std::vector<std::vector<int>> nn_map = octree.kNearestNeighbors(vertices, k);

#pragma omp parallel
  {
#pragma omp for
    for (int i = 0; i < num_vertices; i++) {
      std::vector<int> ni = nn_map[i];
      std::set<int> local_set(ni.begin(), ni.end());
      _adj_list[i] = local_set;
    };
  }
}
