#ifndef UTILS_HPP
#define UTILS_HPP

#include "Octree.hpp"
#include <random>

template <int min = 0, int max = 100>
Octree *rand_tree(std::mt19937 &gen, int num_points, int max_depth) {
  std::uniform_int_distribution<> dis((double)min, (double)max);

  std::vector<std::array<double, 3>> points(num_points);
  for (int i = 0; i < num_points; i++) {
    std::generate(points[i].begin(), points[i].end(),
                  [&]() { return dis(gen); });
  }
  return new Octree(points, max_depth);
};

double distance(std::array<double, 3> a, Node *node);
double distance(std::array<double, 3> a, std::array<double, 3> b);

#endif
