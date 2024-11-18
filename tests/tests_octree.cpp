#include "../src/Octree.hpp"
#include <array>
#include <gtest/gtest.h>
#include <vector>

TEST(OctreeConstruct, Basic) {

  //           .----.----. (1,1,1)
  //          /|   /|   / |
  // (0,1,0) .----.----.__.
  //         |/|  | |  | /|
  //         .----.----. -.
  //         |/   |/   | /
  //         .----.----.
  //        0           (1,0,0)

  std::vector<std::array<double, 3>> points = {
      {0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {1, 1, 1}};
  Octree tree(points);
};
