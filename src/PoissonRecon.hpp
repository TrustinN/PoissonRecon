#ifndef POISSON_RECON_HPP
#define POISSON_RECON_HPP

#include "HRefine.hpp"
#include "PPolynomialXd.hpp"
#include "pOctree.hpp"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <array>
#include <vector>

class PoissonRecon {
public:
  PoissonRecon(const std::vector<std::array<double, 3>> &points,
               const std::vector<std::array<double, 3>> &normals,
               const std::vector<std::array<double, 3>> &inward_normals,
               int depth = 8);

  void run();
  void run_full_solve();
  void write();

  pOctree octree() { return _octree; };
  auto field_centers() const { return _field_centers; };
  auto field_normals() const { return _field_normals; };
  void computeVectorField();
  double evaluateDivergence(const std::array<double, 3> &point);
  double getIsoValue(const std::array<double, 3> &p) const;

private:
  int _depth;
  std::vector<std::array<double, 3>> _points;
  std::vector<std::vector<std::array<double, 3>>> _centers;
  std::vector<std::array<double, 3>> _normals;
  std::vector<std::array<double, 3>> _inward_normals;

  std::vector<std::array<double, 3>> _field_normals;
  std::vector<std::array<double, 3>> _field_centers;

  pOctree _octree;

  std::vector<std::vector<double>> _coeff;
  std::vector<std::vector<ScalarField<2>>> _scalar_fields;
  PPolynomialXD<2, 3> _indicator_function;
};

#endif
