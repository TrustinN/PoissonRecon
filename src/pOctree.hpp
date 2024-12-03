#ifndef P_OCTREE_HPP
#define P_OCTREE_HPP

#include "Octree.hpp"
#include <set>

class pOctree : public Octree {
public:
  pOctree() : Octree() {};
  pOctree(std::vector<std::array<double, 3>> points, int depth = 8);
  void AssignVecField(std::vector<std::array<double, 3>> normals);
  // For each node, we find all 26 neighboring nodes of the same depth within
  // distance d
  std::vector<Node *> Neighbors(Node *node);

  // Compute the inner product of the divergence of our vector field with the
  // basis function on our node
  double ExtractInnerProduct(Node *node);

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

#endif
