#ifndef EMST_HPP
#define EMST_HPP

#include "Octree.hpp"
#include <set>
#include <vector>

class Emst {
  std::vector<std::set<int>> _adj_list;
  int _num_edges;

public:
  Emst(const std::vector<std::array<double, 3>> vertices, Octree<Node> &octree);
  int num_edges() { return _num_edges; };
  std::vector<std::set<int>> adj_list() { return _adj_list; }
};

#endif
