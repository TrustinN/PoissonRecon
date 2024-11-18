#ifndef IO_HPP
#define IO_HPP

#include "Octree.hpp"
#include <array>

std::ostream &operator<<(std::ostream &ofs, const std::array<double, 3> &a);
std::ostream &operator<<(std::ostream &ofs, const Node &n);
std::ostream &operator<<(std::ostream &ofs, const Octree &o);

#endif
