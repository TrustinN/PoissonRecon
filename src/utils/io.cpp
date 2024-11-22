#include "io.hpp"
#include "../Node.hpp"
#include "../Octree.hpp"
#include <array>
#include <iostream>

std::ostream &operator<<(std::ostream &ofs, const std::array<double, 3> &a) {
  ofs << "[" << a[0] << ", " << a[1] << ", " << a[2] << "]";
  return ofs;
};
