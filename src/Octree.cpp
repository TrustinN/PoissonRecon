#include "Octree.hpp"

oNode::oNode(std::vector<id_point> points, std::array<double, 3> center,
             double width, bool is_leaf, int depth)
    : Node(points, center, width, is_leaf, depth) {
  // info is union type, choose if we want to store points or children
  if (is_leaf) {
    new (&info.points) std::vector<id_point>(points);
  } else {
    new (&info.children) std::array<oNode *, 8>{nullptr};
  }
};
