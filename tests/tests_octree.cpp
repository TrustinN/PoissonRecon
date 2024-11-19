#include "../src/Octree.hpp"
#include "../src/utils.hpp"
#include <array>
#include <gtest/gtest.h>
#include <random>
#include <vector>

std::random_device rd;
std::mt19937 gen(rd());

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
  Octree *tree = rand_tree(gen, 100, 2);
  ASSERT_EQ(tree->size(), 100);
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
  Octree *tree = rand_tree(gen, max_points, 8);

  std::array<double, 3> query = {1.0, 2.0, 3.0};

  std::vector<int> result = tree->kNearestNeighbors(query);

  // brute force find answer
  std::array<double, 3> answer;
  double min_dist = std::numeric_limits<double>::infinity();
  for (auto point : tree->points()) {
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

  double nearest_dist = pow(query[0] - tree->points()[n_id][0], 2) +
                        pow(query[1] - tree->points()[n_id][1], 2) +
                        pow(query[2] - tree->points()[n_id][2], 2);
  EXPECT_LE(nearest_dist, answer_dist);
}
