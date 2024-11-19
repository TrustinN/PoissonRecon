#ifndef UTILS_HPP
#define UTILS_HPP

#include "Octree.hpp"

double distance(std::array<double, 3> a, Node *node);
double distance(std::array<double, 3> a, std::array<double, 3> b);

std::vector<int> rand_ints(int min, int max, int num);
std::vector<std::array<double, 3>> rand_points(double min, double max,
                                               int num_points);

#endif
