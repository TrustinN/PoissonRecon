#ifndef BASIS_HPP
#define BASIS_HPP

#include <array>

// 1-D basis function
struct basisF1 {
  static double operator()(double p);
};

// assigns a weight over the range -1 < x < 1 for
// x, y, z
struct basisF {
  double operator()(const std::array<double, 3> &p);
};

#endif
