#include "Normal.hpp"
#include "PoissonRecon.hpp"
#include "basis.hpp"
#include "utils/io.hpp"
#include "utils/linalg.hpp"
#include "utils/sampling.hpp"

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
  // std::vector<std::array<double, 3>> vertices = sample_sphere(1000, 3);
  // NormalApproximations na(vertices);
  // PoissonRecon poisson(vertices, na.normals(), na.inward_normals());
  // std::cout << poisson.v() << std::endl;
  divVField dV;
  std::cout << dV.int_field_x << std::endl;
  std::cout << dV.int_field_y << std::endl;
  std::cout << dV.int_field_z << std::endl;
  std::cout << dV.wc << std::endl;
  std::cout << dV.dw << std::endl;
  //
  laplaceField lf;
  std::cout << lf.int_field_x << std::endl;
  std::cout << lf.int_field_y << std::endl;
  std::cout << lf.int_field_z << std::endl;
  std::cout << lf.wc << std::endl;
  std::cout << lf.dw << std::endl;

  return 0;
}
