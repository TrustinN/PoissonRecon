#ifndef UTILS_HPP
#define UTILS_HPP

#include "Octree.hpp"
#include <random>

Octree *rand_tree(std::mt19937 &gen, int num_points, int max_depth);

#endif
