#include "Node.hpp"
#include "utils/io.hpp"
#include <iostream>

Node::Node(std::vector<id_point> points, std::array<double, 3> center,
           double width, bool is_leaf, int depth)
    : center(center), width(width), depth(depth), is_leaf(is_leaf),
      num_points(points.size()) {};

Node::~Node() {
  if (is_leaf) {
    info.points.~vector();
  } else {
    info.children.~array();
  }
}

std::ostream &operator<<(std::ostream &ofs, const id_point &a) {
  return ofs << std::get<1>(a);
};

std::ostream &operator<<(std::ostream &ofs, const Node &n) {
  std::string indent(2 * n.depth, ' ');
  std::string label = n.is_leaf ? "Leaf" : "Branch";
  ofs << indent << label << " {" << std::endl;
  ofs << indent << "  " << "center: " << n.center << "," << std::endl;
  ofs << indent << "  " << "width: " << n.width << "," << std::endl;
  if (n.is_leaf) {
    ofs << indent << "  " << n.info.points << std::endl;
  } else {
    ofs << indent << "  " << "children: " << n.info.children << std::endl;
  }
  ofs << indent << "}";
  return ofs;
};

std::ostream &operator<<(std::ostream &ofs, const std::array<Node *, 8> &a) {
  if (a[0] == nullptr) {
    ofs << "[null";
  } else {
    ofs << std::endl;
    ofs << *a[0];
  }
  for (int i = 1; i < 8; i++) {
    if (a[i] == nullptr) {
      ofs << ", null";
    } else {
      ofs << std::endl << *a[i];
    }
  };
  if (a[0] == nullptr) {
    ofs << "]";
  }
  return ofs;
};
