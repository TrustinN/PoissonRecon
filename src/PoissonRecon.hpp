#ifndef POISSON_RECON_HPP
#define POISSON_RECON_HPP

#include "Eigen/Dense"
#include "Eigen/Sparse"
#include "basis.hpp"
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
  Eigen::VectorXd v() { return _v; };

private:
  int _depth;
  std::vector<std::array<double, 3>> _points;
  std::vector<std::array<double, 3>> _normals;
  std::vector<std::array<double, 3>> _inward_normals;
  pOctree _octree;
  divergenceField _divergence_field;
  laplacianField _laplacian_field;

  // The system to be solved is _L * _x = _v
  Eigen::SparseMatrix<double, Eigen::ColMajor> _L;
  Eigen::VectorXd _x;
  Eigen::VectorXd _v;
};

#endif
