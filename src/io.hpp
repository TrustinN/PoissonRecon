#ifndef IO_HPP
#define IO_HPP

#include "Octree.hpp"
#include <array>
#include <iostream>

std::ostream &operator<<(std::ostream &ofs, const std::array<double, 3> &a);
std::ostream &operator<<(std::ostream &ofs, const id_point &a);
template <typename T>
std::ostream &operator<<(std::ostream &ofs, const std::vector<T> &v) {
  ofs << "<";
  if (!v.empty()) {
    ofs << v[0];
    for (std::size_t i = 1; i < v.size(); i++) {
      ofs << ", " << v[i];
    };
  }
  ofs << ">";
  return ofs;
}

std::ostream &operator<<(std::ostream &ofs, const Node &n);
std::ostream &operator<<(std::ostream &ofs, const Octree &o);

#endif
