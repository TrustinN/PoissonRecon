#include "sampling.hpp"
#include <random>

std::vector<std::array<double, 3>> sample_sphere(int n, double r) {

  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<double> dis(0.0, r);
  std::uniform_real_distribution<double> dis_scale(0.0, 1.0);

  std::vector<std::array<double, 3>> points(n);
  for (int i = 0; i < n; i++) {
    double z = 2 * r * dis_scale(gen) - r;
    double theta = 2 * M_PI * dis_scale(gen);

    double a = std::sqrt(std::pow(r, 2) - std::pow(z, 2));
    double x = std::cos(theta) * a;
    double y = std::sin(theta) * a;

    points[i] = {x, y, z};
  };
  return points;
};
