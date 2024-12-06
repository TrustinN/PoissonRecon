#ifndef OCTREE_HPP
#define OCTREE_HPP

#include "Node.hpp"
#include "utils/linalg.hpp"
#include <array>
#include <iostream>
#include <queue>
#include <vector>

struct oNode : public Node {
  oNode(std::vector<id_point> points, std::array<double, 3> center,
        double width, bool is_leaf, int depth);
  NodeInfo<oNode> info;
};

template <typename TreeNode> class Octree {
public:
  Octree() : _size(0), _root(nullptr) {};
  Octree(std::vector<std::array<double, 3>> points, int max_depth = 8,
         int min_depth = -1);

  std::vector<int> kNearestNeighbors(const std::array<double, 3> query,
                                     int k = 1) const;
  void Delete(std::array<double, 3> p);

  int size() const { return _size; };

  TreeNode *root() const { return _root; };
  std::vector<int> unused() const { return unused_ids; };
  std::vector<std::array<double, 3>> points() const { return _points; };
  std::vector<TreeNode *> child_nodes() const { return _child_nodes; };

private:
  int _size;
  int _max_depth;
  int _min_depth;
  int _curr_depth;
  TreeNode *_root;
  std::vector<std::array<double, 3>> _points;
  std::vector<int>
      unused_ids; // When we delete, deleted point's ids will go here. This way,
                  // we don't need to resize our _points, when we insert a new
                  // point, we can just replace the value in that position
  std::vector<TreeNode *> _child_nodes;

  TreeNode *build(std::vector<id_point> points, std::array<double, 3> center,
                  double width, int depth);

  int Delete(TreeNode *node, std::array<double, 3> p);
};

// -------------------------------------------------------------------------------------------------//
// Template Instantiation
// -------------------------------------------------------------------------------------------------//
//

template <typename TreeNode>
TreeNode *Octree<TreeNode>::build(std::vector<id_point> points,
                                  std::array<double, 3> center, double width,
                                  int depth) {

  bool is_leaf;
  if (_min_depth == -1) {
    is_leaf = depth == _max_depth || points.size() <= 1;
  } else {
    is_leaf = (_min_depth <= depth && points.size() <= 1) ||
              depth == _max_depth || points.size() == 0;
  };
  TreeNode *ret_node = new TreeNode(points, center, width, is_leaf, depth);
  _curr_depth = std::max(depth, _curr_depth);

  if (!is_leaf) {

    // Figure out which points go in which subdivision
    std::vector<std::vector<id_point>> sub_d(8);
    for (const auto &p : points) {
      std::array<double, 3> cur_p = std::get<1>(p);
      std::array<double, 3> diff = {cur_p[0] - center[0], cur_p[1] - center[1],
                                    cur_p[2] - center[2]};
      diff[0] = (diff[0] > 0) ? 1 : 0;
      diff[1] = (diff[1] > 0) ? 2 : 0;
      diff[2] = (diff[2] > 0) ? 4 : 0;

      int idx = diff[0] + diff[1] + diff[2];

      sub_d[idx].push_back(p);
    };

    // Make pointers to subdivision nodes
    std::array<TreeNode *, 8> &n_child = ret_node->info.children;
    double n_width = width / 2.0;

    for (int i = 0; i < 8; i++) {

      std::array<double, 3> n_center = center;
      n_center[0] += (i % 2 == 0) ? -n_width : n_width;
      n_center[1] += ((i >> 1) % 2 == 0) ? -n_width : n_width;
      n_center[2] += ((i >> 2) % 2 == 0) ? -n_width : n_width;

      n_child[i] = build(sub_d[i], n_center, n_width, depth + 1);
    };
  } else {
    _child_nodes.push_back(ret_node);
  };

  return ret_node;
};

template <typename TreeNode>
Octree<TreeNode>::Octree(std::vector<std::array<double, 3>> points,
                         int max_depth, int min_depth)
    : _size(points.size()), _points(points), _max_depth(max_depth),
      _min_depth(min_depth) {

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

    std::vector<id_point> id_points(points.size());
    for (int i = 0; i < points.size(); i++) {
      id_points[i] = std::make_tuple(i, points[i]);
    }

    this->_root = Octree<TreeNode>::build(id_points, center, width / 2, 1);

  } else {
    this->_root = nullptr;
  };
};

// priority queue data stores either node or point
// based on distance to query
template <typename TreeNode> struct pqData {
  double priority;
  bool is_point;
  int id;

  union data {
    TreeNode *node;
    std::array<double, 3> pt;
    data(TreeNode *node) : node(node) {};
    data(const std::array<double, 3> pt) : pt(pt) {};
    ~data() {};
  } data;

  pqData(double p, TreeNode *node)
      : priority(p), data(node), is_point(false) {};
  pqData(double p, const std::array<double, 3> pt, int id)
      : priority(p), data(pt), is_point(true), id(id) {};

  friend bool operator<(const pqData &lhs, const pqData &rhs) {
    return lhs.priority < rhs.priority;
  }

  friend bool operator>(const pqData &lhs, const pqData &rhs) {
    return lhs.priority > rhs.priority;
  }
};

// helper function computes distance between
template <typename TreeNode>
double distance(const std::array<double, 3> &a, const TreeNode *node) {
  std::array<double, 3> center = node->center;
  std::array<double, 3> diff = {std::abs(a[0] - center[0]),
                                std::abs(a[1] - center[1]),
                                std::abs(a[2] - center[2])};

  // If distance along an axis is < node width
  // only need to compute remaining distance along other axis
  double width = node->width;
  diff[0] = (diff[0] < width) ? 0 : std::pow(diff[0] - width, 2);
  diff[1] = (diff[1] < width) ? 0 : std::pow(diff[1] - width, 2);
  diff[2] = (diff[2] < width) ? 0 : std::pow(diff[2] - width, 2);

  return diff[0] + diff[1] + diff[2];
};

template <typename TreeNode>
std::vector<int>
Octree<TreeNode>::kNearestNeighbors(const std::array<double, 3> query,
                                    int k) const {
  std::vector<int> ret;

  std::priority_queue<pqData<TreeNode>, std::vector<pqData<TreeNode>>,
                      std::greater<pqData<TreeNode>>>
      pq;
  pq.push(pqData(0, this->_root));

  while (ret.size() < k && !pq.empty()) {
    pqData item = pq.top();
    pq.pop();

    if (item.is_point) {
      // pq item is a point, we add id as a nearest neighbor
      ret.push_back(item.id);

    } else {
      // add each child of node to pq by distance to query

      // two cases: if node is a leaf or if node is not a leaf
      TreeNode *node = item.data.node;

      if (node->is_leaf) {
        // we add the points of the leaf
        for (const id_point &p : node->info.points) {
          int id = std::get<0>(p);
          std::array<double, 3> cur_p = std::get<1>(p);
          pq.push(pqData<TreeNode>(distance(query, cur_p), cur_p, id));
        }
      } else {

        // we add children node to pq
        for (const auto &child : node->info.children) {
          if (child != nullptr) {
            pq.push(pqData<TreeNode>(distance(query, child), child));
          }
        }
      }
    }
  };

  return ret;
};

template <typename TreeNode>
int Octree<TreeNode>::Delete(TreeNode *node, std::array<double, 3> p) {
  int del_id = -1;
  if (node->is_leaf) {

    std::vector<id_point> &pts = node->info.points;

    for (auto itr = pts.begin(); itr != pts.end(); itr++) {
      id_point id_p = *itr;
      if (p == std::get<1>(id_p)) {
        del_id = std::get<0>(id_p);
        pts.erase(itr);
        node->num_points -= 1;
        break;
      };
    };

  } else {
    // find out which child node the point is in
    std::array<double, 3> &center = node->center;
    std::array<double, 3> diff = {p[0] - center[0], p[1] - center[1],
                                  p[2] - center[2]};

    diff[0] = (diff[0] > 0) ? 1 : 0;
    diff[1] = (diff[1] > 0) ? 2 : 0;
    diff[2] = (diff[2] > 0) ? 4 : 0;

    int idx = diff[0] + diff[1] + diff[2];
    TreeNode *next_node = node->info.children[idx];
    del_id = Delete(next_node, p);

    if (del_id != -1) {
      node->num_points -= 1;

      // now check if there is only one point left. If so, we can turn this
      // node to a leaf node
      if (node->num_points == 1) {
        id_point l_pt;
        node->is_leaf = true;
        for (int i = 0; i < 8; i++) {
          if (node->info.children[i]->num_points == 0) {
            // clean up
            delete node->info.children[i];
            node->info.children[i] = nullptr;
          } else {
            l_pt = node->info.children[i]->info.points[0];
          }
        }
        // change union data type
        node->info.children.~array();
        new (&node->info.points) std::vector<id_point>{l_pt};

      } else if (node->num_points == 0) {
        node->is_leaf = true;

        for (int i = 0; i < 8; i++) {
          // clean up
          delete node->info.children[i];
          node->info.children[i] = nullptr;
        }
      }
    };
  };

  return del_id;
};

template <typename TreeNode>
void Octree<TreeNode>::Delete(std::array<double, 3> p) {
  int del_id = Delete(_root, p);
  if (del_id != -1) {
    _size -= 1;
    unused_ids.push_back(del_id);
  };
};

#endif
