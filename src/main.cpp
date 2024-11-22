#include "Normal.hpp"
#include "Octree.hpp"
#include "utils/io.hpp"
#include "utils/sampling.hpp"
#include <iostream>

int main() {

  // std::vector<std::array<double, 3>> vertices = sample_sphere(100, 3);
  // NormalApproximations normals(vertices);
  // std::cout << normals.normals() << std::endl;
  // std::vector<std::array<double, 3>> vertices = {
  //     {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0},
  //     {2.0, 0.0, 0.0}, {0.0, 2.0, 0.0}, {0.0, 0.0, 2.0}, {1.0, 1.0, 0.0},
  //     {1.0, 0.0, 1.0}, {0.0, 1.0, 1.0}, {1.0, 1.0, 1.0}, {3.0, 5.0, 3.0},
  //     {2.0, 4.0, 2.0}, {1.0, 3.0, 1.0}};
  std::vector<std::array<double, 3>> vertices = rand_points(1, 100, 1000);

  Octree octree(vertices, 3, -1);
  std::cout << octree << std::endl;

  for (int i = 0; i < vertices.size(); i++) {
    octree.Delete(vertices[i]);
  }
  std::cout << octree << std::endl;

  return 0;
}
