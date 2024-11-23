#include "Octree.hpp"
#include "utils/linalg.hpp"
#include <array>
#include <cassert>
#include <iostream>
#include <queue>
#include <vector>

// -------------------------------------------------------------------------------------------------//
// Helpers
// -------------------------------------------------------------------------------------------------//

// maps the current point to which child node it belongs
int node_index_map(Node *node, const std::array<double, 3> &p) {
  std::array<double, 3> diff = p - node->center;
  diff[0] = (diff[0] > 0) ? 1 : 0;
  diff[1] = (diff[1] > 0) ? 2 : 0;
  diff[2] = (diff[2] > 0) ? 4 : 0;

  return diff[0] + diff[1] + diff[2];
}

std::vector<std::array<double, 3>> split_centers(Node *node) {

  double n_width = node->width / 2.0;
  std::array<double, 3> &center = node->center;

  std::vector<std::array<double, 3>> n_centers;

  for (int i = 0; i < 8; i++) {
    std::array<double, 3> n_center = center;
    n_center[0] += (i % 2 == 0) ? -n_width : n_width;
    n_center[1] += ((i >> 1) % 2 == 0) ? -n_width : n_width;
    n_center[2] += ((i >> 2) % 2 == 0) ? -n_width : n_width;
    n_centers.push_back(n_center);
  };

  return n_centers;
};

std::array<std::vector<id_point>, 8>
partition_points(Node *node, const std::vector<id_point> &points) {
  // set where each point goes
  std::array<std::vector<id_point>, 8> point_partition;
  for (const auto &p : points) {
    int idx = node_index_map(node, std::get<1>(p));
    point_partition[idx].push_back(p);
  };
  return point_partition;
};

// -------------------------------------------------------------------------------------------------//
// Implementation
// -------------------------------------------------------------------------------------------------//

Node::Node(std::vector<id_point> points, std::array<double, 3> center,
           double width, bool is_leaf, int depth)
    : center(center), width(width), depth(depth), is_leaf(is_leaf),
      num_points(points.size()) {
  // info is union type, choose if we want to store points or children
  if (is_leaf) {
    new (&info.points) std::vector<id_point>(points);
  } else {
    new (&info.children) std::array<Node *, 8>{nullptr};
  }
};

Node::~Node() {
  if (is_leaf) {
    info.points.~vector();
  } else {
    info.children.~array();
  }
};

void Node::Insert(const std::vector<id_point> &points) {
  assert(this->is_leaf);
  this->info.points.insert(this->info.points.end(), points.begin(),
                           points.end());
}

void Node::subdivide() {
  assert(this->is_leaf);
  this->is_leaf = false;

  // get current points
  std::vector<id_point> points = std::move(this->info.points);

  // switch union type
  this->info.points.~vector();
  this->info.children = std::array<Node *, 8>{nullptr};

  std::array<std::vector<id_point>, 8> point_partition =
      partition_points(this, points);
  std::vector<std::array<double, 3>> n_centers = split_centers(this);

  // populate children as Node pointers
  for (int i = 0; i < 8; i++) {
    this->info.children[i] = new Node(point_partition[i], n_centers[i],
                                      width / 2, true, this->depth + 1);
  };
};

// -------------------------------------------------------------------------------------------------//
// Tree Implementation
// -------------------------------------------------------------------------------------------------//

Node *Octree::build(std::vector<id_point> points, std::array<double, 3> center,
                    double width, int depth) {

  bool is_leaf;
  if (_min_depth == -1) {
    is_leaf = depth == _max_depth || points.size() <= 1;
  } else {
    is_leaf = (_min_depth <= depth && points.size() == 1) ||
              depth == _max_depth || points.size() == 0;
  };
  Node *ret_node = new Node(points, center, width, is_leaf, depth);

  if (!is_leaf) {

    // Figure out which points go in which subdivision
    std::array<std::vector<id_point>, 8> sub_d =
        partition_points(ret_node, points);

    // Make pointers to subdivision nodes
    std::array<Node *, 8> &n_child = ret_node->info.children;
    std::vector<std::array<double, 3>> n_centers = split_centers(ret_node);

    for (int i = 0; i < 8; i++) {
      n_child[i] = build(sub_d[i], n_centers[i], width / 2.0, depth + 1);
    };
  } else {
    _base_nodes.push_back(ret_node);
  };

  return ret_node;
};

Octree::Octree(std::vector<std::array<double, 3>> points, int max_depth,
               int min_depth)
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

    // expand the bounds by 4 * smallest node width
    double pad = width / std::pow(2, max_depth - 2);
    this->_root = Octree::build(id_points, center, width / 2 + pad, 1);

  } else {
    this->_root = nullptr;
  };
};

// priority queue data stores either node or point
// based on distance to query
struct pqData {
  double priority;
  bool is_point;
  int id;

  union data {
    Node *node;
    std::array<double, 3> pt;
    data(Node *node) : node(node) {};
    data(const std::array<double, 3> pt) : pt(pt) {};
    ~data() {};
  } data;

  pqData(double p, Node *node) : priority(p), data(node), is_point(false) {};
  pqData(double p, const std::array<double, 3> pt, int id)
      : priority(p), data(pt), is_point(true), id(id) {};

  friend bool operator<(const pqData &lhs, const pqData &rhs) {
    return lhs.priority < rhs.priority;
  }

  friend bool operator>(const pqData &lhs, const pqData &rhs) {
    return lhs.priority > rhs.priority;
  }
};

std::vector<int> Octree::kNearestNeighbors(const std::array<double, 3> query,
                                           int k) const {
  std::vector<int> ret;

  std::priority_queue<pqData, std::vector<pqData>, std::greater<pqData>> pq;
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
      Node *node = item.data.node;

      if (node->is_leaf) {
        // we add the points of the leaf
        for (const id_point &p : node->info.points) {
          int id = std::get<0>(p);
          std::array<double, 3> cur_p = std::get<1>(p);
          pq.push(pqData(distance(query, cur_p), cur_p, id));
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

int Octree::Delete(Node *node, std::array<double, 3> p) {
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
    int idx = node_index_map(node, p);
    Node *next_node = node->info.children[idx];
    del_id = Delete(next_node, p);

    if (del_id != -1) {
      node->num_points -= 1;

      if (node->num_points == 0) {
        node->is_leaf = true;
        for (int i = 0; i < 8; i++) {
          // clean up
          delete node->info.children[i];
          node->info.children[i] = nullptr;
        }

        // change union data type
        node->info.children.~array();
      }
    };
  };

  return del_id;
};

void Octree::Delete(std::array<double, 3> p) {
  int del_id = Delete(_root, p);
  if (del_id != -1) {
    _size -= 1;
    unused_ids.push_back(del_id);
  };
};
