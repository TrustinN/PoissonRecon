#ifndef POISSON_RECON_HPP
#define POISSON_RECON_HPP

#include "pOctree.hpp"
#include <array>
#include <vector>

class PoissonRecon {
public:
  PoissonRecon(const std::vector<std::array<double, 3>> &points,
               const std::vector<std::array<double, 3>> &normals,
               const std::vector<std::array<double, 3>> &inward_normals,
               int depth = 8);

  pOctree octree() { return _octree; };

private:
  int _depth;
  std::vector<std::array<double, 3>> _points;
  std::vector<std::array<double, 3>> _normals;
  std::vector<std::array<double, 3>> _inward_normals;
  pOctree _octree;
};

#endif
