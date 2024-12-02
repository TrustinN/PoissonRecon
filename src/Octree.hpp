#ifndef OCTREE_HPP
#define OCTREE_HPP

#include <array>
#include <cassert>
#include <vector>

using id_point = std::pair<int, std::array<double, 3>>;

// -------------------------------------------------------------------------------------------------//
// Node
// -------------------------------------------------------------------------------------------------//

//      6.----.----. 7      \\ index into children by bit map
//      /|   /|  3/ |       \\ xyz where x, y, z are 0 or 1
//    2.----.----.__.
//     |/|4 | |  | /|
//     .----.----. -.5
//     |/   |/   | /
//     .----.----.
//    0          1
struct Node;

union NodeChildren {
  std::vector<id_point> points;
  std::array<Node *, 8> nodes;
  NodeChildren() : points() {};
  ~NodeChildren() {};
};

struct Node {
  Node(std::vector<id_point> points, std::array<double, 3> center, double width,
       bool is_leaf, int depth);
  ~Node();

  void Insert(const std::vector<id_point> &points);
  void subdivide();

  double width;
  std::array<double, 3> center;

  int depth;
  int depth_id = -1;

  bool is_leaf;
  int num_points;
  NodeChildren children;

  std::array<double, 3> normal = {0.0, 0.0, 0.0};
};

// -------------------------------------------------------------------------------------------------//
// Octree
// -------------------------------------------------------------------------------------------------//

class Octree {

protected:
  int _size;
  Node *_root;
  int _max_depth;
  int _min_depth;

  std::vector<int> _deleted_point_ids;
  std::vector<std::vector<int>> _deleted_node_ids;

  std::vector<std::array<double, 3>> _points;
  std::vector<std::vector<Node *>> _nodes;

  Node *build(std::vector<id_point> points, std::array<double, 3> center,
              double width, int depth);

  int Delete(Node *node, std::array<double, 3> p);
  void register_node(Node *node);
  void unregister_node(Node *node);

private:
  std::vector<int> &getUnusedNodeIds(int d) {
    assert(0 <= d && d <= _deleted_node_ids.size());
    return _deleted_node_ids[d];
  };

public:
  Octree() : _size(0), _root(nullptr) {};
  Octree(std::vector<std::array<double, 3>> points, int max_depth = 8,
         int min_depth = -1);

  std::vector<int> kNearestNeighbors(const std::array<double, 3> query,
                                     int k = 1) const;
  template <bool refine = false>
  void Insert(const std::vector<std::array<double, 3>> &points);
  template <bool refine = false>
  void Insert(Node *node, const std::vector<id_point> &points);
  void Delete(std::array<double, 3> p);

  std::vector<Node *> &getNodesAtDepth(int d) {
    assert(0 <= d && d <= _nodes.size());
    return _nodes[d];
  };
  std::vector<Node *> getNodesAtDepth(int d) const {
    assert(0 <= d && d <= _nodes.size());
    return _nodes[d];
  };

  Node *root() const { return _root; };
  int size() const { return _size; };
  std::vector<int> deleted_ids() const { return _deleted_point_ids; };
  std::vector<std::array<double, 3>> points() const { return _points; };
};

// -------------------------------------------------------------------------------------------------//
// Helper Functions
// -------------------------------------------------------------------------------------------------//

std::vector<std::array<double, 3>> split_centers(Node *node);
int node_index_map(Node *node, const std::array<double, 3> &p);
std::array<std::vector<id_point>, 8>
partition_points(Node *node, const std::vector<id_point> &points);

// -------------------------------------------------------------------------------------------------//
// Template Implementations
// -------------------------------------------------------------------------------------------------//

template <bool refine>
void Octree::Insert(Node *node, const std::vector<id_point> &points) {

  if (node->depth_id < 0) {
    this->register_node(node);
  }

  if (points.size() == 0) {
    return;
  } else if (!refine) {
    node->num_points += points.size();
  }

  if (!node->is_leaf || node->depth < _min_depth) {
    if (node->is_leaf) {
      node->subdivide();
      node->is_leaf = false;
    }

    // Figure out which points go in which subdivision
    auto point_partition = partition_points(node, points);

    for (int i = 0; i < 8; i++) {
      this->Insert<refine>(node->children.nodes[i], point_partition[i]);
    }
    return;
  }

  // we are at the target depth and the node is a leaf node
  if (!refine) {
    node->Insert(points);
  }
}

template <bool refine>
void Octree::Insert(const std::vector<std::array<double, 3>> &points) {
  // Assumes that we do not need to expand Octree bounds
  std::vector<id_point> id_points(points.size());
  int start_id = _points.size();
  for (int i = 0; i < points.size(); i++) {
    id_points[i] = {start_id + i, points[i]};
  }
  this->Insert<refine>(_root, id_points);

  if (!refine) {
    // update points
    _points.insert(_points.end(), points.begin(), points.end());
    _size += points.size();
  }
};

#endif
