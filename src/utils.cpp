#include "utils.hpp"
#include "Octree.hpp"
#include <array>

// helper function computes distance between
double distance(std::array<double, 3> a, Node *node) {
  std::array<double, 3> center = node->center;
  std::array<double, 3> diff = {std::abs(a[0] - center[0]),
                                std::abs(a[1] - center[1]),
                                std::abs(a[2] - center[2])};

  // If distance along an axis is < node width
  // only need to compute remaining distance along other axis
  double width = node->width / 2;
  diff[0] = (diff[0] < width) ? 0 : std::pow(diff[0] - width, 2);
  diff[1] = (diff[1] < width) ? 0 : std::pow(diff[1] - width, 2);
  diff[2] = (diff[2] < width) ? 0 : std::pow(diff[2] - width, 2);

  return diff[0] + diff[1] + diff[2];
};

double distance(std::array<double, 3> a, std::array<double, 3> b) {
  return std::pow(a[0] - b[0], 2) + std::pow(a[1] - b[1], 2) +
         std::pow(a[2] - b[2], 2);
};
