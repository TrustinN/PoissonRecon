#ifndef EMST_HPP
#define EMST_HPP

#include "Octree.hpp"
#include <set>
#include <vector>

class Emst {
  std::vector<std::set<int>> adj_list;
  int num_edges;

  Emst(const std::vector<std::array<double, 3>> vertices, Octree &octree);
};

#endif
