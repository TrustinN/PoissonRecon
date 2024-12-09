
#ifndef GRAPHS_HPP
#define GRAPHS_HPP

#include <queue>
#include <set>

// -------------------------------------------------------------------------------------------------//
// GRAPH UTILS
// -------------------------------------------------------------------------------------------------//

template <typename T> using weight = double(const T &a, const T &b);

struct MSTpqData {
  double priority;
  int id_root;
  int id_leaf;

  MSTpqData(){};
  MSTpqData(double p, int id_root, int id_leaf)
      : priority(p), id_root(id_root), id_leaf(id_leaf){};

  friend bool operator>(const MSTpqData &lhs, const MSTpqData &rhs) {
    return lhs.priority > rhs.priority;
  };
};

template <typename T>
std::vector<std::set<int>> get_mst(const std::vector<T> &vertices,
                                   const std::vector<std::set<int>> &adj_list,
                                   weight<T> weight) {
  int edges = 0;
  std::vector<std::set<int>> ret_adj_list(vertices.size());
  std::vector<int> visited(vertices.size(), 0);
  visited[0] = 1;

  std::priority_queue<MSTpqData, std::vector<MSTpqData>,
                      std::greater<MSTpqData>>
      min_pq;

  for (int idx : adj_list[0]) {
    min_pq.push(MSTpqData(weight(vertices[0], vertices[idx]), 0, idx));
  };

  while (true) {
    MSTpqData item = min_pq.top();
    min_pq.pop();
    if (!visited[item.id_leaf]) {
      visited[item.id_leaf] = 1;

      ret_adj_list[item.id_root].insert(item.id_leaf);
      ret_adj_list[item.id_leaf].insert(item.id_root);
      edges++;

      if (edges == vertices.size() - 1) {
        break;
      }

      for (int idx : adj_list[item.id_leaf]) {
        if (!visited[idx]) {
          min_pq.push(MSTpqData(weight(vertices[item.id_leaf], vertices[idx]),
                                item.id_leaf, idx));
        };
      };
    };
  };

  return ret_adj_list;
}

std::vector<std::set<int>> join_graphs(const std::vector<std::set<int>> &g1,
                                       const std::vector<std::set<int>> &g2);
#endif
