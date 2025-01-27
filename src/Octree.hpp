#ifndef OCTREE_HPP
#define OCTREE_HPP

#include "Metrics.hpp"
#include <array>
#include <cassert>
#include <vector>

using id_point = std::pair<int, std::array<double, 3>>;

// -------------------------------------------------------------------------------------------------//
// Node
// -------------------------------------------------------------------------------------------------//

//      6.----.----. 7
//      /|   /|  3/ |
//    2.----.----.__.
//     |/|4 | |  | /|
//     .----.----. -.5
//     |/   |/   | /
//     .----.----.
//    0          1
// index into children by bit map
// xyz where x, y, z are 0 or 1
//

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

  std::vector<int>
  kNearestNeighbors(const std::array<double, 3> &query, int k = 1,
                    const Metric &metric = DefaultMetric()) const;

  std::vector<std::vector<int>>
  kNearestNeighbors(const std::vector<std::array<double, 3>> &queries,
                    int k) const;

  std::vector<Node *> Refine(const std::vector<std::array<double, 3>> &points);
  void Refine(Node *node, const std::vector<std::array<double, 3>> &points,
              std::vector<Node *> &refinement);
  void Delete(std::array<double, 3> p);
  std::vector<int> RadiusSearch(const std::array<double, 3> &center, double r,
                                int depth);
  std::vector<Node *> Neighbors(Node *node);

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
  int node_count() const;
  int max_depth() const { return _max_depth; };
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
std::array<std::vector<std::array<double, 3>>, 8>
partition_points(Node *node, const std::vector<std::array<double, 3>> &points);
Node *seek_node(Node *node, const std::array<double, 3> &p);
Node *seek_node(Node *start, const std::array<double, 3> &p, int depth);
std::vector<std::array<double, 3>> nearest_8(Node *node,
                                             const std::array<double, 3> &p);
std::vector<std::array<double, 3>> nearest_27(Node *node);

// -------------------------------------------------------------------------------------------------//
// Template Implementations
// -------------------------------------------------------------------------------------------------//

double distance(const std::array<double, 3> &a, const Node *node);

std::ostream &operator<<(std::ostream &ofs, const id_point &a);
std::ostream &operator<<(std::ostream &ofs, const Node &n);
std::ostream &operator<<(std::ostream &ofs, Node *n);
std::ostream &operator<<(std::ostream &ofs, const Octree &o);

#endif
