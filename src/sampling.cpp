#include "sampling.hpp"
#include <random>

std::vector<std::array<double, 3>> sample_sphere(int n, double r) {

  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<double> dis(0.0, r);
  std::uniform_real_distribution<double> dis_scale(0.0, 1.0);

  std::vector<std::array<double, 3>> points(n);
  for (int i = 0; i < n; i++) {
    double x = dis(gen);
    double y = std::sqrt(std::pow(r, 2) - std::pow(x, 2));

    double scale = dis_scale(gen);

    x *= scale;
    y *= scale;

    double base = r * scale;
    double height = std::sqrt(std::pow(r, 2) - std::pow(base, 2));

    points[i] = {x, y, height};
  };
  return points;
};
