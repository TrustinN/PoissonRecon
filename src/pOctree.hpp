#ifndef P_OCTREE_HPP
#define P_OCTREE_HPP

#include "Node.hpp"
#include "Octree.hpp"

// 1-D basis function
struct basisF1 {
  static double operator()(double p);
};

// assigns a weight over the range -1 < x < 1 for
// x, y, z
struct basisF {
  double operator()(const std::array<double, 3> &p);
};

struct pNode : public Node {
  pNode(std::vector<id_point> points, std::array<double, 3> center,
        double width, bool is_leaf, int depth);
  NodeInfo<pNode> info;
};

class pOctree : public Octree<pNode> {
public:
  pOctree() : Octree<pNode>() {};
  pOctree(std::vector<std::array<double, 3>> points, int depth = 8)
      : Octree<pNode>(points, depth, depth) {};
};

#endif
