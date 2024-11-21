#include "../src/Octree.hpp"
#include "../src/utils/utils.hpp"
#include <array>
#include <gtest/gtest.h>
#include <numeric>
#include <vector>

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

TEST(OctreeConstruct, Random) {
  std::vector<std::array<double, 3>> points = rand_points(0, 100, 1000);
  Octree tree(points);
  ASSERT_EQ(tree.size(), 1000);
}

TEST(OctreekNN, Basic) {
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
  std::vector<int> nn = tree.kNearestNeighbors({0.2, 0, 0});
  int nn_id = nn[0];

  ASSERT_EQ(nn.size(), 1);
  ASSERT_EQ(nn_id, 0);

  nn = tree.kNearestNeighbors({0.7, 0.7, 0.7});
  nn_id = nn[0];
  ASSERT_EQ(nn.size(), 1);
  ASSERT_EQ(nn_id, 3);
}

TEST(OctreekNN, Random) {

  const int max_points = 10000;
  std::vector<std::array<double, 3>> points = rand_points(0, 100, max_points);
  Octree tree(points);

  std::array<double, 3> query = {1.0, 2.0, 3.0};

  std::vector<int> result = tree.kNearestNeighbors(query);

  // brute force find answer
  std::array<double, 3> answer;
  double min_dist = std::numeric_limits<double>::infinity();
  for (auto point : tree.points()) {
    std::array<double, 3> diff;
    std::transform(point.begin(), point.end(), query.begin(), diff.begin(),
                   std::minus<double>());
    double sD = std::accumulate(
        diff.begin(), diff.end(), 0.0,
        [](double sum, double a) { return sum + std::pow(a, 2); });
    if (sD < min_dist) {
      min_dist = sD;
      answer = point;
    }
  }

  int n_id = result[0];

  double answer_dist = pow(query[0] - answer[0], 2) +
                       pow(query[1] - answer[1], 2) +
                       pow(query[2] - answer[2], 2);

  double nearest_dist = pow(query[0] - tree.points()[n_id][0], 2) +
                        pow(query[1] - tree.points()[n_id][1], 2) +
                        pow(query[2] - tree.points()[n_id][2], 2);
  EXPECT_LE(nearest_dist, answer_dist);
};

TEST(OctreeDelete, NoDelete) {
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
  tree.Delete({1, 2, 3});

  ASSERT_EQ(tree.size(), 4);
  ASSERT_EQ(tree.unused().size(), 0);
};

TEST(OctreeDelete, Basic) {

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
  tree.Delete({0, 0, 0});
  ASSERT_EQ(tree.size(), 3);
  ASSERT_EQ(tree.unused()[0], 0);
  ASSERT_EQ(tree.unused().size(), 1);
};

TEST(OctreeDelete, FullDelete) {

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
  tree.Delete({0, 0, 0});
  tree.Delete({1, 0, 0});
  tree.Delete({0, 1, 0});
  tree.Delete({1, 1, 1});
  ASSERT_EQ(tree.size(), 0);
  ASSERT_EQ(tree.unused().size(), 4);
  ASSERT_EQ(tree.root()->is_leaf, true);
  ASSERT_EQ(tree.root()->info.points.size(), 0);
};

TEST(OctreeDelete, Random) {
  std::vector<std::array<double, 3>> points = rand_points(0, 100, 1000);
  Octree tree(points);
  std::vector<int> indices = rand_ints(0, 200, 1000);

  for (int idx : indices) {
    tree.Delete(tree.points()[idx]);
  }

  std::set<int> expected_del_ids = std::set(indices.begin(), indices.end());
  std::set<int> actual_del_ids =
      std::set(tree.unused().begin(), tree.unused().end());

  ASSERT_EQ(expected_del_ids, actual_del_ids);
}
