#include "Octree.hpp"
#include "pOctree.hpp"
#include "utils/io.hpp"
#include "utils/sampling.hpp"
#include <iostream>

int main() {

  // std::vector<std::array<double, 3>> vertices = {
  //     {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0},
  //     {2.0, 0.0, 0.0}, {0.0, 2.0, 0.0}, {0.0, 0.0, 2.0}, {1.0, 1.0, 0.0},
  //     {1.0, 0.0, 1.0}, {0.0, 1.0, 1.0}, {1.0, 1.0, 1.0}, {3.0, 5.0, 3.0},
  //     {2.0, 4.0, 2.0}, {1.0, 3.0, 1.0}};
  std::vector<std::array<double, 3>> vertices = rand_points(1, 100, 20);

  Octree octree(vertices, 4, 4);
  std::cout << octree << std::endl;

  // pOctree octree2({}, 4);
  // octree2.Insert(vertices);
  // std::cout << octree2 << std::endl;

  pOctree octree2(vertices, 4);
  std::cout << octree2 << std::endl;

  return 0;
}
