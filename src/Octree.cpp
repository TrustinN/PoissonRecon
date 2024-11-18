#include "Octree.hpp"
#include <array>
#include <vector>

// -------------------------------------------------------------------------------------------------//
// Implementation
// -------------------------------------------------------------------------------------------------//

Node::Node(std::vector<std::array<double, 3>> points,
           std::array<double, 3> center, double width, bool is_leaf)
    : center(center), width(width) {
  // info is union type, choose if we want to store points or children
  if (is_leaf) {
    new (&info.points) std::vector<std::array<double, 3>>(points);
  } else {
    new (&info.children) std::array<Node *, 8>{nullptr};
  }
};

Node *Octree::build(std::vector<std::array<double, 3>> points,
                    std::array<double, 3> center, double width, int depth) {
  Node *ret_node = new Node(points, center, width, depth == 0);
  if (points.size() == 0) {
    return ret_node;
  };

  if (depth > 0) {

    // Figure out which points go in which subdivision
    std::vector<std::vector<std::array<double, 3>>> sub_d(8);
    for (const auto &p : points) {

      std::array<double, 3> diff;
      std::transform(center.begin(), center.end(), p.begin(), diff.begin(),
                     std::minus<double>());
      std::transform(diff.begin(), diff.end(), diff.begin(),
                     [](double &d) { return (d > 0) ? 1 : 0; });

      int idx = (int)diff[0] + ((int)diff[1] << 1) + ((int)diff[2] << 2);
      sub_d[idx].push_back(p);
    };

    // Make pointers to subdivision nodes
    std::array<Node *, 8> n_child = ret_node->info.children;
    double n_width = width / 2.0;

    for (int i = 0; i < 8; i++) {

      std::array<double, 3> n_center = center;
      n_center[0] += (i % 2 == 0) ? -n_width : n_width;
      n_center[1] += ((i >> 1) % 2 == 0) ? -n_width : n_width;
      n_center[2] += ((i >> 2) % 2 == 0) ? -n_width : n_width;

      n_child[i] = build(sub_d[i], n_center, n_width, depth - 1);
    };
  };
  return ret_node;
};

Octree::Octree(std::vector<std::array<double, 3>> points, int max_depth) {
  this->_size = points.size();
  // compute center and width of containing tree

  std::array<double, 3> mins{std::numeric_limits<double>::infinity()};
  std::array<double, 3> maxes{-std::numeric_limits<double>::infinity()};

  for (const auto &p : points) {
    std::transform(
        mins.begin(), mins.end(), p.begin(), mins.begin(),
        [](const double &a, const double &b) { return std::min(a, b); });
    std::transform(
        maxes.begin(), maxes.end(), p.begin(), maxes.begin(),
        [](const double &a, const double &b) { return std::max(a, b); });
  };

  std::array<double, 3> center = {(maxes[0] + mins[0]) / 2.0,
                                  (maxes[1] + mins[1]) / 2.0,
                                  (maxes[2] + mins[2]) / 2.0};
  double width = std::max(maxes[0] - mins[0],
                          std::max(maxes[1] - mins[1], maxes[2] - mins[2]));

  this->_root = Octree::build(points, center, width, max_depth);
};
