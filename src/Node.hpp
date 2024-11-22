#ifndef NODE_HPP
#define NODE_HPP

#include <vector>

using id_point = std::pair<int, std::array<double, 3>>;

template <typename TreeNode> union NodeInfo {
  std::vector<id_point> points;
  std::array<TreeNode *, 8> children;
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

  double width;
  std::array<double, 3> center;
  bool is_leaf;
  int depth;
  int num_points;
  NodeInfo<Node> info;
};

std::ostream &operator<<(std::ostream &ofs, const id_point &a);
std::ostream &operator<<(std::ostream &ofs, const Node &n);
std::ostream &operator<<(std::ostream &ofs, const std::array<Node *, 8> &a);

#endif
