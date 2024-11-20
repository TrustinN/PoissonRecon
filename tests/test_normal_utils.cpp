#include "../src/Normal.hpp"
#include "../src/utils.hpp"
#include <array>
#include <gtest/gtest.h>
#include <set>
#include <vector>

TEST(NormalUtilsTP, Basic) {
  std::vector<std::array<double, 3>> vertices = {
      {0, 0, 0}, {0, 0, 1}, {0, 1, 0}};
  std::array<double, 3> normal = get_normal(vertices);
  std::array<double, 3> expected_normal = {1, 0, 0};
  ASSERT_EQ(normal, expected_normal);
}

TEST(NormalUtilsMST, Basic) {
  //               .(0, 2, 2)
  //               |
  //               |
  //               |    <--- Expected MST
  //               |
  // (0, 0, 0)     |
  // .-------------.(0, 0, 2)

  std::vector<std::array<double, 3>> vertices = {
      {0, 0, 0}, {0, 0, 2}, {0, 2, 2}};
  std::vector<std::set<int>> adj_list = {
      {1, 2}, {0, 2}, {0, 1}}; // Fully connected

  std::vector<std::set<int>> mst =
      get_mst<std::array<double, 3>>(vertices, adj_list, distance);
  std::vector<std::set<int>> expected_mst = {{1}, {0, 2}, {1}};
  ASSERT_EQ(mst, expected_mst);
}
