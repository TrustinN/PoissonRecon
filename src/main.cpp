#include "Normal.hpp"
// #include "pOctree.hpp"
#include "Node.hpp"
#include "Octree.hpp"
#include "utils/io.hpp"
#include "utils/sampling.hpp"
#include <iostream>

int main() {

  // std::vector<std::array<double, 3>> vertices = sample_sphere(100, 3);
  // NormalApproximations normals(vertices);
  // std::cout << normals.normals() << std::endl;
  std::vector<std::array<double, 3>> vertices = {
      {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};
  // Octree<Node> octree(vertices);
  // std::cout << octree << std::endl;
  Octree<oNode> octree(vertices);
  // for (const auto &v : octree.child_nodes()) {
  //   std::cout << *v << std::endl;
  // }
  // std::cout << octree << std::endl;

  return 0;
}
