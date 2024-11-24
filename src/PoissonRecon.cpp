#include "PoissonRecon.hpp"
#include <iostream>

PoissonRecon::PoissonRecon(
    const std::vector<std::array<double, 3>> &points,
    const std::vector<std::array<double, 3>> &normals,
    const std::vector<std::array<double, 3>> &inward_normals, int depth)
    : _points(points), _normals(normals), _inward_normals(inward_normals),
      _depth(depth) {

  _octree = pOctree(points, depth);
  _octree.AssignVecField(normals);

  // compute v
  for (Node *node : _octree.field_nodes()) {
    _v.push_back(_octree.ExtractInnerProduct(node));
  }
};
