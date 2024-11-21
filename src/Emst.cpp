#include "Emst.hpp"
#include "utils/graphs.hpp" // Contains MSTpqData
#include "utils/linalg.hpp" // Distance function
#include <queue>

Emst::Emst(const std::vector<std::array<double, 3>> vertices, Octree &octree)
    : _num_edges(0) {

  _adj_list = std::vector<std::set<int>>(vertices.size());

  std::vector<int> visited(vertices.size(), 0);
  visited[0] = 1;

  // stores id and distance to tree
  std::priority_queue<MSTpqData, std::vector<MSTpqData>,
                      std::greater<MSTpqData>>
      min_pq;
  octree.Delete(vertices[0]);
  int id_new = octree.kNearestNeighbors(vertices[0])[0];

  min_pq.push(MSTpqData(0, 0, id_new));
  while (true) {
    MSTpqData item = min_pq.top();
    min_pq.pop();
    if (!visited[item.id_leaf]) {
      visited[item.id_leaf] = 1;

      // add leaf to tree
      _adj_list[item.id_root].insert(item.id_leaf);
      _adj_list[item.id_leaf].insert(item.id_root);
      _num_edges++;

      if (_num_edges == vertices.size() - 1) {
        break;
      }

      // get corresponding vertices
      std::array<double, 3> v_root = vertices[item.id_root];
      std::array<double, 3> v_leaf = vertices[item.id_leaf];

      // delete item from octree so nearest neighbor finds a non mst node
      octree.Delete(v_leaf);

      // push both items and their new nearest neighbors to the tree
      int v_root_id_leaf = octree.kNearestNeighbors(v_root)[0];
      int v_leaf_id_leaf = octree.kNearestNeighbors(v_leaf)[0];

      min_pq.push(MSTpqData(distance(v_root, vertices[v_root_id_leaf]),
                            item.id_root, v_root_id_leaf));
      min_pq.push(MSTpqData(distance(v_leaf, vertices[v_leaf_id_leaf]),
                            item.id_leaf, v_leaf_id_leaf));
    };
  };
};
