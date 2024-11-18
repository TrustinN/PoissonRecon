#include "Octree.hpp"
#include "io.hpp"
#include "utils.hpp"
#include <iostream>
#include <random>

int main() {
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<> dis(0.0, 100.0);

  Octree *tree = rand_tree(gen, 1000, 8);
  std::cout << *tree << std::endl;

  return 0;
}
