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

  std::vector<std::array<double, 3>> field_centers() const {
    return _field_centers;
  };
  std::vector<std::array<double, 3>> field_normals() const {
    return _field_normals;
  };

private:
  std::vector<std::array<double, 3>> _field_centers;
  std::vector<std::array<double, 3>> _field_normals;
};

Node *seek_node(Node *node, const std::array<double, 3> &p);
std::vector<std::array<double, 3>> nearest_8(Node *node,
                                             const std::array<double, 3> &p);
std::vector<std::array<double, 3>> nearest_27(Node *node);

#endif
