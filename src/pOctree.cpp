#include "pOctree.hpp"

double basisF1::operator()(double p) {
  double abs_p = std::abs(p);
  if (abs_p >= 1.5) {
    return 0;
  } else if (abs_p > .5) {
    return std::pow(1.5 - abs_p, 2) / 2;
  } else {
    return .75 - std::pow(p, 2);
  };
}

double basisF::operator()(const std::array<double, 3> &p) {
  return basisF1::operator()(p[0]) * basisF1::operator()(p[1]) *
         basisF1::operator()(p[2]);
};

pNode::pNode(std::vector<id_point> points, std::array<double, 3> center,
             double width, bool is_leaf, int depth)
    : Node(points, center, width, is_leaf, depth) {
  // info is union type, choose if we want to store points or children
  if (is_leaf) {
    new (&info.points) std::vector<id_point>(points);
  } else {
    new (&info.children) std::array<pNode *, 8>{nullptr};
  }
};
