#include "linalg.hpp"

// -------------------------------------------------------------------------------------------------//
// LINEAR ALGEBRA
// -------------------------------------------------------------------------------------------------//

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
