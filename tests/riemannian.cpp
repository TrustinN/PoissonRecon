// #include "../src/RiemannianGraph.hpp"
// #include "../src/utils/sampling.hpp"
// #include <gtest/gtest.h>
//
// TEST(RiemannianGraphConstruct, Parallel) {
//   const int max_points = 10000;
//   const int nn = 15;
//   std::vector<std::array<double, 3>> points =
//       rand_points(-100, 100, max_points);
//   Octree tree(points);
//   RiemannianGraph rG1(points, tree);
//   RiemannianGraph rG2(points, tree, true);
//
//   std::vector<std::set<int>> adj1 = rG1.adj_list();
//   std::vector<std::set<int>> adj2 = rG2.adj_list();
//   for (int i = 0; i < points.size(); i++) {
//     ASSERT_EQ(adj1[i], adj2[i]);
//   }
// }
