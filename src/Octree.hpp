#ifndef OCTREE_HPP
#define OCTREE_HPP

#include <array>
#include <iostream>
#include <vector>

using id_point = std::pair<int, std::array<double, 3>>;

struct Node;

union NodeInfo {
  std::vector<id_point> points;
  std::array<Node *, 8> children;
  NodeInfo() : points() {};
  ~NodeInfo() {};
};

struct Node {
  //      6.----.----. 7      \\ index into children by bit map
  //      /|   /|  3/ |       \\ xyz where x, y, z are 0 or 1
  //    2.----.----.__.
  //     |/|4 | |  | /|
  //     .----.----. -.5
  //     |/   |/   | /
  //     .----.----.
  //    0          1
  Node(std::vector<id_point> points, std::array<double, 3> center, double width,
       bool is_leaf, int depth);
  ~Node();

  void Insert(const std::vector<id_point> &points);
  void subdivide();

  double width;
  std::array<double, 3> center;
  bool is_leaf;
  int depth;
  int num_points;
  NodeInfo info;

  // not native to octree implementation
  // too lazy to extend from base class
  std::array<double, 3> normal = {0.0, 0.0, 0.0};
  int id = -1;
};

class Octree {
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

  int size() const { return _size; };
  Node *root() const { return _root; };
  std::vector<std::array<double, 3>> points() const { return _points; };
  std::vector<int> deleted_ids() const { return _deleted_ids; };
  std::vector<Node *> leaf_nodes() const { return _leaf_nodes; };
  std::vector<int> unused_leaves() const { return _unused_leaves; };

protected:
  int _size;
  int _max_depth;
  int _min_depth;
  std::vector<Node *> _leaf_nodes;
  Node *_root;
  std::vector<std::array<double, 3>> _points;
  std::vector<int> _deleted_ids; // When we delete, deleted point's ids will go
                                 // here. This way, we don't need to resize our
                                 // _points, when we insert a new point, we can
                                 // just replace the value in that position
  std::vector<int>
      _unused_leaves; // When we delete, deleted point's ids will go
                      // here. This way, we don't need to resize our
                      // _points, when we insert a new point, we can
                      // just replace the value in that position

  Node *build(std::vector<id_point> points, std::array<double, 3> center,
              double width, int depth);

  int Delete(Node *node, std::array<double, 3> p);
  void insert_as_leaf(Node *node);
  void remove_leaf(Node *node);
};

std::vector<std::array<double, 3>> split_centers(Node *node);
int node_index_map(Node *node, const std::array<double, 3> &p);
std::array<std::vector<id_point>, 8>
partition_points(Node *node, const std::vector<id_point> &points);

template <bool refine>
void Octree::Insert(Node *node, const std::vector<id_point> &points) {

  if (node->is_leaf) {
    // add as leaf node
    if (node->id < 0) {
      this->insert_as_leaf(node);
    }
  }

  if (points.size() > 0) {
    if (!refine) {
      node->num_points += points.size();
    }

    if (node->is_leaf) {
      if (node->depth >= _min_depth) {

        if (!refine) {
          node->Insert(points);
        }
      } else {
        // we must increase the depth
        std::array<std::vector<id_point>, 8> point_partition =
            partition_points(node, points);

        node->subdivide();
        this->remove_leaf(node);

        for (int i = 0; i < 8; i++) {
          this->Insert<refine>(node->info.children[i], point_partition[i]);
        };
      };
    } else {

      // Figure out which points go in which subdivision
      std::array<std::vector<id_point>, 8> point_partition =
          partition_points(node, points);

      for (int i = 0; i < 8; i++) {
        this->Insert<refine>(node->info.children[i], point_partition[i]);
      }
    }
  }
};

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
