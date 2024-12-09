#include "Normal.hpp"
#include "PoissonRecon.hpp"
#include "utils/sampling.hpp"

int main(int argc, char **argv) {
  std::string prop = "ball";
  for (int i = 1; i < argc; i++) {
    std::string input = std::string(argv[i]);
    if (input == "--box") {
      prop = "box";
    }
  }
  std::vector<std::array<double, 3>> vertices;
  if (prop == "ball") {
    vertices = sample_sphere(5000, 20);
  } else {
    vertices = sample_box(5000, 15, 20, 10);
  }
  NormalApproximations na(vertices);
  PoissonRecon poisson(vertices, na.normals(), na.inward_normals(), 6);
  poisson.run();
  poisson.write();

  return 0;
}
