#include "Normal.hpp"
#include "PoissonRecon.hpp"
#include "utils/sampling.hpp"

int main(int argc, char **argv) {
  std::string prop = "ball";
  int num_samples = -1;
  for (int i = 1; i < argc; i++) {
    std::string input = std::string(argv[i]);
    if (input == "--ball") {
      prop = "ball";
      num_samples = std::stoi(argv[i + 1]);
    }
    if (input == "--box") {
      prop = "box";
      num_samples = std::stoi(argv[i + 1]);
    }
  }
  if (num_samples < -1) {
    std::cout << "Unspecified object to sample" << std::endl;
    return 0;
  }

  std::vector<std::array<double, 3>> vertices;
  if (prop == "ball") {
    vertices = sample_sphere(num_samples, 20);
  } else {
    vertices = sample_box(num_samples, 15, 20, 10);
  }
  NormalApproximations na(vertices);
  PoissonRecon poisson(vertices, na.normals(), na.inward_normals(), 6);
  poisson.run();
  poisson.write();

  return 0;
}
