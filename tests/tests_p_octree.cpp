#include "../src/pOctree.hpp"
#include "../src/utils/io.hpp"
#include "../src/utils/sampling.hpp"
#include <array>
#include <gtest/gtest.h>
#include <numeric>
#include <vector>

TEST(pOctreeConstruct, Empty) {
  pOctree *tree = new pOctree(std::vector<std::array<double, 3>>(), 3);
  ASSERT_EQ(tree->size(), 0);
};

TEST(pOctreeConstruct, Default) {
  pOctree tree;
  ASSERT_EQ(tree.size(), 0);
};

TEST(pOctreeConstruct, Basic) {

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
  pOctree tree(points);
  ASSERT_EQ(tree.size(), 4);
};

TEST(pOctreeConstruct, Random) {
  std::vector<std::array<double, 3>> points = rand_points(0, 100, 1000);
  pOctree tree(points);
  std::vector<pNode *> child_nodes = tree.child_nodes();

  for (pNode *child : child_nodes) {
    for (int i = 0; i < child->info.points.size(); i++) {
      id_point id_p = child->info.points[i];
      ASSERT_EQ(tree.points()[std::get<0>(id_p)], std::get<1>(id_p));
    }
  }

  ASSERT_EQ(tree.size(), 1000);
}

TEST(pOctreekNN, Basic) {
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
  pOctree tree(points);
  std::vector<int> nn = tree.kNearestNeighbors({0.2, 0, 0});
  int nn_id = nn[0];

  ASSERT_EQ(nn.size(), 1);
  ASSERT_EQ(nn_id, 0);

  nn = tree.kNearestNeighbors({0.7, 0.7, 0.7});
  nn_id = nn[0];
  ASSERT_EQ(nn.size(), 1);
  ASSERT_EQ(nn_id, 3);
}

TEST(pOctreekNN, Random) {

  const int max_points = 10000;
  std::vector<std::array<double, 3>> points = rand_points(0, 100, max_points);
  pOctree tree(points);

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

TEST(pOctreeDelete, NoDelete) {
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
  pOctree tree(points);
  tree.Delete({1, 2, 3});

  ASSERT_EQ(tree.size(), 4);
  ASSERT_EQ(tree.unused().size(), 0);
};

TEST(pOctreeDelete, Basic) {

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
  pOctree tree(points);
  tree.Delete({0, 0, 0});
  ASSERT_EQ(tree.size(), 3);
  ASSERT_EQ(tree.unused()[0], 0);
  ASSERT_EQ(tree.unused().size(), 1);
};

TEST(pOctreeDelete, FullDelete) {

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
  pOctree tree(points);

  std::vector<pNode *> child_nodes = tree.child_nodes();
  for (pNode *child : child_nodes) {
    for (int i = 0; i < child->info.points.size(); i++) {
      id_point id_p = child->info.points[i];
      std::cout << std::get<1>(id_p) << std::endl;
    }
  }
  tree.Delete({0, 0, 0});
  tree.Delete({1, 0, 0});
  tree.Delete({0, 1, 0});
  tree.Delete({1, 1, 1});
  // ASSERT_EQ(tree.size(), 0);
  // ASSERT_EQ(tree.unused().size(), 4);
  ASSERT_EQ(tree.root()->is_leaf, true);
  ASSERT_EQ(tree.root()->info.points.size(), 0);
};

TEST(pOctreeDelete, Random) {
  std::vector<std::array<double, 3>> points = rand_points(0, 100, 1000);
  pOctree tree(points);
  std::vector<int> indices = rand_ints(0, 200, 1000);

  for (int idx : indices) {
    tree.Delete(tree.points()[idx]);
  }

  std::set<int> expected_del_ids = std::set(indices.begin(), indices.end());
  std::set<int> actual_del_ids =
      std::set(tree.unused().begin(), tree.unused().end());

  ASSERT_EQ(expected_del_ids, actual_del_ids);
}
