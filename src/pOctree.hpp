#ifndef P_OCTREE_HPP
#define P_OCTREE_HPP

#include "Octree.hpp"

// 1-D basis function
struct basisF1 {
  static double operator()(double p);
};

// assigns a weight over the range -1 < x < 1 for
// x, y, z
struct basisF {
  double operator()(const std::array<double, 3> &p);
};

class pOctree : Octree {
public:
  pOctree() : Octree() {};
  pOctree(std::vector<std::array<double, 3>> points, int depth = 8)
      : Octree(points, depth, depth) {};

private:
};

#endif
