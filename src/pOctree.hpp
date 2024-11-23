#ifndef P_OCTREE_HPP
#define P_OCTREE_HPP

#include "Octree.hpp"
#include <set>

class pOctree : public Octree {
public:
  pOctree() : Octree() {};
  pOctree(std::vector<std::array<double, 3>> points, int depth = 8);
  void AssignVecField(std::vector<std::array<double, 3>> normals);
  std::set<std::array<double, 3>> field_loc() { return _field_centers; };
  int center_count() { return _field_centers.size(); };
  std::vector<Node *> field_nodes() { return _field_nodes; };

private:
  std::set<std::array<double, 3>> _field_centers;
  std::vector<Node *> _field_nodes;
};

#endif
