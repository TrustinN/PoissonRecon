#include "Octree.hpp"
#include <array>
#include <numeric>
#include <queue>
#include <vector>

// -------------------------------------------------------------------------------------------------//
// Implementation
// -------------------------------------------------------------------------------------------------//

Node::Node(std::vector<std::array<double, 3>> points,
           std::array<double, 3> center, double width, bool is_leaf, int depth)
    : center(center), width(width), depth(depth), is_leaf(is_leaf) {
  // info is union type, choose if we want to store points or children
  if (is_leaf) {
    new (&info.points) std::vector<std::array<double, 3>>(points);
  } else {
    new (&info.children) std::array<Node *, 8>{nullptr};
  }
};

Node *Octree::build(std::vector<std::array<double, 3>> points,
                    std::array<double, 3> center, double width, int depth,
                    int max_depth) {

  bool is_leaf = depth == max_depth || points.size() <= 1;
  Node *ret_node = new Node(points, center, width, is_leaf, depth);

  if (!is_leaf) {

    // Figure out which points go in which subdivision
    std::vector<std::vector<std::array<double, 3>>> sub_d(8);
    for (const auto &p : points) {

      std::array<double, 3> diff;
      std::transform(p.begin(), p.end(), center.begin(), diff.begin(),
                     std::minus<double>());
      std::transform(diff.begin(), diff.end(), diff.begin(),
                     [](double &d) { return (d > 0) ? 1 : 0; });

      int idx = (int)diff[0] + ((int)diff[1] << 1) + ((int)diff[2] << 2);
      sub_d[idx].push_back(p);
    };

    // Make pointers to subdivision nodes
    std::array<Node *, 8> &n_child = ret_node->info.children;
    double n_width = width / 2.0;
    double c_off = n_width / 2.0;

    for (int i = 0; i < 8; i++) {

      std::array<double, 3> n_center = center;
      n_center[0] += (i % 2 == 0) ? -c_off : c_off;
      n_center[1] += ((i >> 1) % 2 == 0) ? -c_off : c_off;
      n_center[2] += ((i >> 2) % 2 == 0) ? -c_off : c_off;

      n_child[i] = build(sub_d[i], n_center, n_width, depth + 1, max_depth);
    };
  };

  return ret_node;
};

Octree::Octree(std::vector<std::array<double, 3>> points, int max_depth)
    : _size(points.size()), _points(points) {

  if (points.size() > 0) {

    // compute center and width of containing tree
    std::array<double, 3> mins{std::numeric_limits<double>::infinity()};
    std::array<double, 3> maxes{-std::numeric_limits<double>::infinity()};

    for (const auto &p : points) {
      std::transform(
          mins.begin(), mins.end(), p.begin(), mins.begin(),
          [](const double &a, const double &b) { return std::min(a, b); });
      std::transform(
          maxes.begin(), maxes.end(), p.begin(), maxes.begin(),
          [](const double &a, const double &b) { return std::max(a, b); });
    };

    std::array<double, 3> center = {(maxes[0] + mins[0]) / 2.0,
                                    (maxes[1] + mins[1]) / 2.0,
                                    (maxes[2] + mins[2]) / 2.0};
    double width = std::max(maxes[0] - mins[0],
                            std::max(maxes[1] - mins[1], maxes[2] - mins[2]));
    this->_root = Octree::build(points, center, width, 1, max_depth);

  } else {
    this->_root = nullptr;
  };
};

// priority queue data stores either node or point
// based on distance to query
struct pqData {
  double priority;
  bool is_point;

  union data {
    Node *node;
    std::array<double, 3> pt;
    data(Node *node) : node(node) {};
    data(const std::array<double, 3> pt) : pt(pt) {};
    ~data() {};
  } data;

  pqData(double p, Node *node) : priority(p), data(node), is_point(false) {};
  pqData(double p, const std::array<double, 3> pt)
      : priority(p), data(pt), is_point(true) {};

  friend bool operator<(const pqData &lhs, const pqData &rhs) {
    return lhs.priority < rhs.priority;
  }

  friend bool operator>(const pqData &lhs, const pqData &rhs) {
    return lhs.priority > rhs.priority;
  }
};

// helper function computes distance between
double distance(std::array<double, 3> a, Node *node) {
  std::array<double, 3> center = node->center;
  std::array<double, 3> diff = {std::abs(a[0] - center[0]),
                                std::abs(a[1] - center[1]),
                                std::abs(a[2] - center[2])};

  // If distance along an axis is < node width
  // only need to compute remaining distance along other axis
  double width = node->width / 2;
  diff[0] = (diff[0] < width) ? 0 : std::pow(diff[0] - width, 2);
  diff[1] = (diff[1] < width) ? 0 : std::pow(diff[1] - width, 2);
  diff[2] = (diff[2] < width) ? 0 : std::pow(diff[2] - width, 2);

  return diff[0] + diff[1] + diff[2];
};

double distance(std::array<double, 3> a, std::array<double, 3> b) {
  return std::pow(a[0] - b[0], 2) + std::pow(a[1] - b[1], 2) +
         std::pow(a[2] - b[2], 2);
};

std::vector<std::array<double, 3>>
Octree::kNearestNeighbors(std::array<double, 3> query, int k) {
  std::vector<std::array<double, 3>> ret;

  std::priority_queue<pqData, std::vector<pqData>, std::greater<pqData>> pq;
  pq.push(pqData(0, this->_root));

  while (ret.size() < k && !pq.empty()) {
    pqData item = pq.top();
    pq.pop();

    if (item.is_point) {
      // pq item is a point, we add this as a nearest neighbor
      ret.push_back(item.data.pt);

    } else {
      // add each child of node to pq by distance to query

      // two cases: if node is a leaf or if node is not a leaf
      Node *node = item.data.node;

      if (node->is_leaf) {
        // we add the points of the leaf
        for (const auto &p : node->info.points) {
          pq.push(pqData(distance(query, p), p));
        }
      } else {

        // we add children node to pq
        for (const auto &child : node->info.children) {
          if (child != nullptr) {
            pq.push(pqData(distance(query, child), child));
          }
        }
      }
    }
  };

  return ret;
};
