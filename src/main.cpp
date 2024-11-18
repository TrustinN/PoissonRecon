#include "Octree.hpp"
#include "io.hpp"
#include "utils.hpp"
#include <iostream>
#include <random>

int main() {
  std::random_device rd;
  std::mt19937 gen(rd());

  Octree *tree = rand_tree(gen, 4, 8);
  std::cout << *tree << std::endl;

  return 0;
}
