#include <array>
#include <vector>

// -------------------------------------------------------------------------------------------------//
// Classes
// -------------------------------------------------------------------------------------------------//

struct Node;

union NodeInfo {
  std::vector<std::array<double, 3>> points;
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
  Node(std::vector<std::array<double, 3>> points, std::array<double, 3> center,
       double width, bool is_leaf);

  double width;
  std::array<double, 3> center;
  bool is_leaf;
  NodeInfo info;
};

class Octree {
public:
  Octree(std::vector<std::array<double, 3>> points, int max_depth = 8);

private:
  int _size;
  Node *_root;
  Node *build(std::vector<std::array<double, 3>> points,
              std::array<double, 3> center, double width, int depth);
};
