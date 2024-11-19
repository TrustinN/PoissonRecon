#include "Octree.hpp"
#include "io.hpp"
#include "utils.hpp"
#include <iostream>
#include <random>

int main() {
  std::random_device rd;
  std::mt19937 gen(rd());

  Octree *tree = rand_tree(gen, 10000, 8);
  // std::cout << *tree << std::endl;

  std::vector<std::array<double, 3>> nn =
      tree->kNearestNeighbors({1.0, 2.0, 3.0}, 3);
  std::cout << nn << std::endl;
  std::cout << nn[0] << std::endl;
  std::cout << nn[1] << std::endl;
  std::cout << nn[2] << std::endl;

  return 0;
}
