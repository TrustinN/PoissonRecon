#ifndef RIEMANNIAN_GRAPH_HPP
#define RIEMANNIAN_GRAPH_HPP

#include "Octree.hpp"
#include <array>
#include <set>
#include <vector>

class RiemannianGraph {
  double _max_edge;
  std::vector<std::set<int>> _adj_list;

public:
  RiemannianGraph(const std::vector<std::array<double, 3>> &vertices,
                  const Octree<oNode> &octree, int k = 15);

  double max_edge() { return _max_edge; };
  std::vector<std::set<int>> adj_list() { return _adj_list; };
};

#endif
