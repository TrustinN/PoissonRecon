#include "graphs.hpp"

// -------------------------------------------------------------------------------------------------//
// GRAPH UTILS
// -------------------------------------------------------------------------------------------------//

std::vector<std::set<int>> join_graphs(const std::vector<std::set<int>> &g1,
                                       const std::vector<std::set<int>> &g2) {
  std::vector<std::set<int>> ret_adj_list(g1.size());
  for (int i = 0; i < g1.size(); i++) {
    std::set_union(g1[i].begin(), g1[i].end(), g2[i].begin(), g2[i].end(),
                   std::inserter(ret_adj_list[i], ret_adj_list[i].begin()));
  }
  return ret_adj_list;
};
