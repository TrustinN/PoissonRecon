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
  std::vector<std::array<double, 3>> vertices = sample_sphere(1500, 3);
  NormalApproximations na(vertices);
  PoissonRecon poisson(vertices, na.normals(), na.inward_normals(), 6);

  // divVField dV;
  // std::string indent = std::string(40, '-');
  // std::cout << indent << "X" << indent << std::endl;
  // std::cout << dV.int_field_x << std::endl;
  // std::cout << dV._infl_x << std::endl;
  // std::cout << indent << "Y" << indent << std::endl;
  // std::cout << dV.int_field_y << std::endl;
  // std::cout << dV._infl_y << std::endl;
  // std::cout << indent << "Z" << indent << std::endl;
  // std::cout << dV.int_field_z << std::endl;
  // std::cout << dV._infl_z << std::endl;
  // std::cout << indent << "Integrals" << indent << std::endl;
  // std::cout << dV.wc << std::endl;
  // std::cout << dV.dw << std::endl;
  // double totl = 0;
  // for (double t : dV.int_field_x) {
  //   totl += t;
  // };
  // std::cout << "sum: " << totl << std::endl;
  //
  // laplaceField lp;
  // std::cout << indent << "X" << indent << std::endl;
  // std::cout << lp.int_field_x << std::endl;
  // std::cout << lp._infl_x << std::endl;
  // std::cout << indent << "Y" << indent << std::endl;
  // std::cout << lp.int_field_y << std::endl;
  // std::cout << lp._infl_y << std::endl;
  // std::cout << indent << "Z" << indent << std::endl;
  // std::cout << lp.int_field_z << std::endl;
  // std::cout << lp._infl_z << std::endl;
  // std::cout << indent << "Integrals" << indent << std::endl;
  // std::cout << lp.wc << std::endl;
  // std::cout << lp.dw << std::endl;
  //
  // totl = 0;
  // for (double t : lp.int_field_x) {
  //   totl += t;
  // };
  // std::cout << "sum: " << totl << std::endl;

  return 0;
}
