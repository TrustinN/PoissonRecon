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

// helper function computes distance between
double distance(const std::array<double, 3> &a, const Node *node) {
  std::array<double, 3> center = node->center;
  std::array<double, 3> diff = {std::abs(a[0] - center[0]),
                                std::abs(a[1] - center[1]),
                                std::abs(a[2] - center[2])};

  // If distance along an axis is < node width
  // only need to compute remaining distance along other axis
  double width = node->width;
  diff[0] = (diff[0] < width) ? 0 : std::pow(diff[0] - width, 2);
  diff[1] = (diff[1] < width) ? 0 : std::pow(diff[1] - width, 2);
  diff[2] = (diff[2] < width) ? 0 : std::pow(diff[2] - width, 2);

  return diff[0] + diff[1] + diff[2];
};

double distance(const std::array<double, 3> &a,
                const std::array<double, 3> &b) {
  return std::pow(a[0] - b[0], 2) + std::pow(a[1] - b[1], 2) +
         std::pow(a[2] - b[2], 2);
};
