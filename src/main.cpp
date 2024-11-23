#include "Normal.hpp"
#include "PoissonRecon.hpp"
#include "utils/sampling.hpp"
#include <iostream>

int main() {

  std::vector<std::array<double, 3>> points = sample_sphere(100, 3);

  std::cout << "Computing normals" << std::endl;
  NormalApproximations na(points);
  std::cout << "Done!" << std::endl;

  std::cout << "Reconstructing surface" << std::endl;
  PoissonRecon poisson(points, na.normals(), na.inward_normals());
  std::cout << "Done!" << std::endl;

  return 0;
}
