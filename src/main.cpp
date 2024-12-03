#include "Normal.hpp"
#include "PoissonRecon.hpp"
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
  std::vector<std::array<double, 3>> vertices = sample_sphere(10000, 1);
  // std::vector<std::array<double, 3>> vertices = sample_box(5000, 1.5, 2, 1);
  NormalApproximations na(vertices);
  PoissonRecon poisson(vertices, na.normals(), na.inward_normals(), 6);
  poisson.run();
  poisson.write();

  // divergenceField df;
  // laplacianField lf;
  //
  // std::cout << df.x_field << std::endl;
  // std::cout << df.y_field << std::endl;
  // std::cout << df.z_field << std::endl;
  // std::cout << lf.x_field << std::endl;
  // std::cout << lf.y_field << std::endl;
  // std::cout << lf.z_field << std::endl;

  // std::string indent = std::string(40, '-');
  // std::array<double, 3> c1 = {1, 1, 1};
  // double width = .5;
  // std::array<double, 3> c2{0};
  // for (int i = 0; i < 27; i++) {
  //   std::array<int, 3> axis = {i % 3, (i / 3) % 3, (i / 9) % 3};
  //   c2[0] = c1[0] + width * (axis[0] - 1);
  //   c2[1] = c1[1] + width * (axis[1] - 1);
  //   c2[2] = c1[2] + width * (axis[2] - 1);
  //
  //   std::cout << indent << i << indent << std::endl;
  //   // std::cout << projection(lf, c1, c2) << std::endl;
  //   // std::cout << projection(lf, c2, c1) << std::endl;
  // };

  return 0;
}
