#include "Octree.hpp"
#include "utils.hpp"
#include <array>
#include <queue>
#include <vector>

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

Node *Octree::build(std::vector<id_point> points, std::array<double, 3> center,
                    double width, int depth, int max_depth) {

  bool is_leaf = depth == max_depth || points.size() <= 1;
  Node *ret_node = new Node(points, center, width, is_leaf, depth);

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

    std::vector<id_point> id_points(points.size());
    for (int i = 0; i < points.size(); i++) {
      id_points[i] = std::make_tuple(i, points[i]);
    }

    this->_root = Octree::build(id_points, center, width, 1, max_depth);

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
    std::array<double, 3> &center = node->center;
    std::array<double, 3> diff = {p[0] - center[0], p[1] - center[1],
                                  p[2] - center[2]};

    diff[0] = (diff[0] > 0 ? 1 : 0);
    diff[1] = (diff[1] > 0 ? 2 : 0);
    diff[2] = (diff[2] > 0 ? 4 : 0);

    int idx = diff[0] + diff[1] + diff[2];
    Node *next_node = node->info.children[idx];
    if (next_node != nullptr) {
      del_id = Delete(next_node, p);

      if (del_id != -1) {
        // Check if child node is empty, if so, delete it
        node->num_points -= 1;
        if (next_node->is_leaf) {
          if (next_node->num_points == 0) {
            // destroy node
            delete next_node;
            node->info.children[idx] = nullptr;
          }
        };

        // now check if there is only one point left. If so, we can turn this
        // node to a leaf node
        if (node->num_points == 1) {
          node->is_leaf = true;
          id_point l_pt;
          for (int i = 0; i < 8; i++) {
            if (node->info.children[i] != nullptr) {
              l_pt = node->info.children[i]->info.points[0];
            }
          }

          // change union data type
          node->info.children.~array();
          new (&node->info.points) std::vector<id_point>{l_pt};
        }
      };
    }; // if node is nullptr, there is nothing to delete
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
