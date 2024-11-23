#ifndef P_OCTREE_HPP
#define P_OCTREE_HPP

#include "Octree.hpp"

class pOctree : public Octree {
public:
  pOctree() : Octree() {};
  pOctree(std::vector<std::array<double, 3>> points, int depth = 8);
};

#endif
