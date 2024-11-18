#include "utils.hpp"
#include "Octree.hpp"
#include <array>
#include <random>

Octree *rand_tree(std::mt19937 &gen, int num_points, int max_depth) {
  std::uniform_int_distribution<> dis(0, 100);

  std::vector<std::array<double, 3>> points(num_points);
  for (int i = 0; i < num_points; i++) {
    std::generate(points[i].begin(), points[i].end(),
                  [&]() { return dis(gen); });
  }
  return new Octree(points, max_depth);
};
