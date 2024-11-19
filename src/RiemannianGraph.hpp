#ifndef RIEMANNIAN_GRAPH_HPP
#define RIEMANNIAN_GRAPH_HPP

#include "Octree.hpp"
#include <array>
#include <set>
#include <vector>

class RiemannianGraph {
  double max_edge;
  std::vector<std::set<int>> adj_list;

  RiemannianGraph(const std::vector<std::array<double, 3>> &vertices,
                  const Octree &octree, int k = 15);
};

#endif
