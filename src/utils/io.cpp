#include "io.hpp"
#include "../Octree.hpp"
#include <array>
#include <iostream>

std::ostream &operator<<(std::ostream &ofs, const std::array<double, 3> &a) {
  ofs << "[" << a[0] << ", " << a[1] << ", " << a[2] << "]";
  return ofs;
};

std::ostream &operator<<(std::ostream &ofs, const id_point &a) {
  return ofs << std::get<1>(a);
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

std::ostream &operator<<(std::ostream &ofs, const Node &n) {
  std::string indent(2 * n.depth, ' ');
  ofs << indent << "Node {" << std::endl;
  ofs << indent << "  " << "center: " << n.center << "," << std::endl;
  ofs << indent << "  " << "width: " << n.width << "," << std::endl;
  if (n.is_leaf) {
    ofs << indent << "  " << n.info.points << std::endl;
  } else {
    ofs << indent << "  " << "children:" << n.info.children << std::endl;
  }
  ofs << indent << "}";
  return ofs;
};

std::ostream &operator<<(std::ostream &ofs, const Octree<Node> &o) {
  ofs << "Octree: [" << std::endl;
  ofs << *o.root() << std::endl;
  ofs << "]";
  return ofs;
}
