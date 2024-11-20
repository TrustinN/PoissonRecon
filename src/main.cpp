#include "Normal.hpp"
#include "io.hpp"
#include "sampling.hpp"
#include <iostream>

int main() {

  std::vector<std::array<double, 3>> vertices = sample_sphere(100, 3);
  NormalApproximations normals(vertices);
  std::cout << normals.normals() << std::endl;

  return 0;
}
