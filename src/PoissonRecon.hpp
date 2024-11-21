#ifndef POISSON_RECON_HPP
#define POISSON_RECON_HPP

#include "pOctree.hpp"
#include <array>
#include <vector>

class PoissonRecon {
public:
  PoissonRecon();

private:
  std::vector<std::array<double, 3>> _points;
  std::vector<std::array<double, 3>> _inward_normals;
  pOctree _octree;

  // for every node, assign to it a vector defining the normal
  // centered at _center
  std::vector<std::array<double, 3>> _v_field;
};

#endif
