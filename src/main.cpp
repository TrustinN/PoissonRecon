#include "Octree.hpp"
#include "io.hpp"
#include "utils.hpp"
#include <iostream>

int main() {
  Octree tree(rand_points(0.0, 100.0, 4), 8);

  std::vector<int> nn = tree.kNearestNeighbors({1.0, 2.0, 3.0}, 3);
  std::cout << nn << std::endl;
  std::cout << nn[0] << std::endl;
  std::cout << nn[1] << std::endl;
  std::cout << nn[2] << std::endl;

  return 0;
}
