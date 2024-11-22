#include "linalg.hpp"

// -------------------------------------------------------------------------------------------------//
// LINEAR ALGEBRA
// -------------------------------------------------------------------------------------------------//

std::array<double, 3> operator*(const std::array<double, 3> &a,
                                const double &b) {
  return std::array<double, 3>{a[0] * b, a[1] * b, a[2] * b};
};

std::array<double, 3> operator*(const double &a,
                                const std::array<double, 3> &b) {
  return b * a;
};

std::array<double, 3> operator+(const std::array<double, 3> &a,
                                const double &b) {
  return std::array<double, 3>{a[0] + b, a[1] + b, a[2] + b};
};

std::array<double, 3> operator+(const double &a,
                                const std::array<double, 3> &b) {
  return b + a;
};

std::array<double, 3> operator+(const std::array<double, 3> &a,
                                const std::array<double, 3> &b) {
  return std::array<double, 3>{a[0] + b[0], a[1] + b[1], a[2] + b[2]};
};

std::array<double, 3> operator-(const std::array<double, 3> &a,
                                const std::array<double, 3> &b) {
  return std::array<double, 3>{a[0] - b[0], a[1] - b[1], a[2] - b[2]};
};

double dot(const std::array<double, 3> &a, const std::array<double, 3> &b) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

double distance(const std::array<double, 3> &a,
                const std::array<double, 3> &b) {
  return std::pow(a[0] - b[0], 2) + std::pow(a[1] - b[1], 2) +
         std::pow(a[2] - b[2], 2);
};
