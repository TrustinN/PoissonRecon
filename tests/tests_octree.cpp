#include "../src/Octree.hpp"
#include "../src/utils.hpp"
#include <array>
#include <gtest/gtest.h>
#include <random>
#include <vector>

// std::random_device rd;
// std::mt19937 gen(rd());
// std::uniform_real_distribution<> dis(0.0, 100.0);

TEST(OctreeConstruct, Empty) {
  Octree *tree = new Octree(std::vector<std::array<double, 3>>(), 3);
  ASSERT_EQ(tree->size(), 0);
};

TEST(OctreeConstruct, Default) {
  Octree tree;
  ASSERT_EQ(tree.size(), 0);
};

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
  ASSERT_EQ(tree.size(), 4);
};

// TEST(OctreeConstruct, Random) {
//   Octree *tree = ::rand_tree(gen, 100, 2);
//   ASSERT_EQ(tree->size(), 100);
// }
