#include "Normal.hpp"
#include "PoissonRecon.hpp"
#include "utils/io.hpp"
#include "utils/sampling.hpp"
#include <iostream>

std::ostream &operator<<(std::ostream &os, const std::array<double, 27> &arr) {
  std::string indent = std::string(3, ' ');
  os << "[" << std::endl << indent << arr[0];

  for (int i = 1; i < 27; i++) {
    os << ", ";
    if (i % 3 == 0) {
      os << std::endl;
    }
    os << indent << arr[i];
  };
  return os << std::endl << "]";
};

int main() {
  std::vector<std::array<double, 3>> vertices = sample_box(5000, 15, 20, 10);
  NormalApproximations na(vertices);
  PoissonRecon poisson(vertices, na.normals(), na.inward_normals(), 6);
  poisson.run();
  poisson.write();

  return 0;
}
