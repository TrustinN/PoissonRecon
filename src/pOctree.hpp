#ifndef P_OCTREE_HPP
#define P_OCTREE_HPP

#include "Octree.hpp"

class pOctree : public Octree {
public:
  pOctree() : Octree() {};
  pOctree(std::vector<std::array<double, 3>> points, int depth = 8);
  void AssignVecField(std::vector<std::array<double, 3>> normals);
  // For each node, we find all 26 neighboring nodes of the same depth within
  // distance d
  std::vector<Node *> Neighbors(Node *node);
  std::vector<int> RadiusSearch(const std::array<double, 3> &center, int r);
};

Node *seek_node(Node *node, const std::array<double, 3> &p);
Node *seek_node(Node *start, const std::array<double, 3> &p, int depth);
std::vector<std::array<double, 3>> nearest_8(Node *node,
                                             const std::array<double, 3> &p);
std::vector<std::array<double, 3>> nearest_27(Node *node);

#endif
