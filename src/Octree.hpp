#ifndef OCTREE_HPP
#define OCTREE_HPP

#include <array>
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

  double width;
  std::array<double, 3> center;
  bool is_leaf;
  int depth;
  NodeInfo info;
};

class Octree {
public:
  Octree() : _size(0), _root(nullptr) {};
  Octree(std::vector<std::array<double, 3>> points, int max_depth = 8);

  std::vector<int> kNearestNeighbors(const std::array<double, 3> query,
                                     int k = 1) const;

  bool Delete(Node *node, std::array<double, 3> p);

  int size() const { return _size; }
  Node *root() const { return _root; }
  std::vector<std::array<double, 3>> points() const { return _points; }

private:
  int _size;
  Node *_root;
  std::vector<std::array<double, 3>> _points;

  Node *build(std::vector<id_point> points, std::array<double, 3> center,
              double width, int depth, int max_depth);
};

#endif
