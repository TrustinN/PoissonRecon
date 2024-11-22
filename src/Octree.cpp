#include "Octree.hpp"
#include <array>
#include <vector>

// -------------------------------------------------------------------------------------------------//
// Implementation
// -------------------------------------------------------------------------------------------------//

Node::Node(std::vector<id_point> points, std::array<double, 3> center,
           double width, bool is_leaf, int depth)
    : center(center), width(width), depth(depth), is_leaf(is_leaf),
      num_points(points.size()) {
  // info is union type, choose if we want to store points or children
  if (is_leaf) {
    new (&info.points) std::vector<id_point>(points);
  } else {
    new (&info.children) std::array<Node *, 8>{nullptr};
  }
};

Node::~Node() {
  if (is_leaf) {
    info.points.~vector();
  } else {
    info.children.~array();
  }
};
