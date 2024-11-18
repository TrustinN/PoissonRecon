#include "io.hpp"
#include "Octree.hpp"
#include <array>
#include <iostream>

std::ostream &operator<<(std::ostream &ofs, const std::array<double, 3> &a) {
  ofs << "[" << a[0] << ", " << a[1] << ", " << a[2] << "]";
  return ofs;
};

std::ostream &operator<<(std::ostream &ofs, const std::array<Node *, 8> &a) {
  ofs << "[";
  if (a[0] == nullptr) {
    ofs << "null";
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
  ofs << "]";
  return ofs;
};

template <typename T>
std::ostream &operator<<(std::ostream &ofs, const std::vector<T> &v) {
  ofs << "<";
  if (!v.empty()) {
    ofs << v[0];
    for (int i = 1; i < v.size(); i++) {
      ofs << ", " << v[i];
    };
  }
  ofs << ">";
  return ofs;
}

std::ostream &operator<<(std::ostream &ofs, const Node &n) {
  std::string indent = "  ";
  for (int i = 0; i < n.depth; i++) {
    ofs << indent;
  }
  ofs << "Node: [" << std::endl;
  for (int i = 0; i < n.depth + 1; i++) {
    ofs << indent;
  }
  ofs << "center: " << n.center << "," << std::endl;
  for (int i = 0; i < n.depth + 1; i++) {
    ofs << indent;
  }
  ofs << "width: " << n.width << "," << std::endl;
  for (int i = 0; i < n.depth + 1; i++) {
    ofs << indent;
  }
  if (n.is_leaf) {
    ofs << n.info.points << std::endl;
  } else {
    ofs << "children:" << n.info.children << std::endl;
  }
  for (int i = 0; i < n.depth; i++) {
    ofs << indent;
  }
  ofs << "]";
  return ofs;
};

std::ostream &operator<<(std::ostream &ofs, const Octree &o) {
  ofs << "Octree: [" << std::endl;
  ofs << *o.root() << std::endl;
  ofs << "]";
  return ofs;
}
