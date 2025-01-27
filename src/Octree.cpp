#include "Octree.hpp"
#include "utils/io.hpp"
#include "utils/linalg.hpp"
#include <array>
#include <cassert>
#include <iostream>
#include <omp.h>
#include <queue>
#include <vector>

constexpr static double EPSILON = 1e-8;

// -------------------------------------------------------------------------------------------------//
// Helpers
// -------------------------------------------------------------------------------------------------//

// helper function computes distance between
double distance(const std::array<double, 3> &a, const Node *node) {
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

std::array<std::vector<std::array<double, 3>>, 8>
partition_points(Node *node, const std::vector<std::array<double, 3>> &points) {
  // set where each point goes
  std::array<std::vector<std::array<double, 3>>, 8> point_partition;
  for (const auto &p : points) {
    int idx = node_index_map(node, p);
    point_partition[idx].push_back(p);
  };
  return point_partition;
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

Node *seek_node(Node *node, const std::array<double, 3> &p) {
  Node *r_node = node;
  while (!r_node->is_leaf) {
    int idx = node_index_map(r_node, p);
    r_node = r_node->children.nodes[idx];
  };
  assert(node->center == p);
  return r_node;
}

Node *seek_node(Node *start, const std::array<double, 3> &p, int depth) {
  Node *r_node = start;
  while (true) {
    if (r_node == nullptr) {
      break;
    }
    if (r_node->depth == depth) {
      if (distance(r_node->center, p) > EPSILON) {
        r_node = nullptr;
      }
      return r_node;
    }
    int idx = node_index_map(r_node, p);
    r_node = r_node->children.nodes[idx];
  };
  return r_node;
}

std::vector<std::array<double, 3>> nearest_8(Node *node,
                                             const std::array<double, 3> &p) {

  std::array<double, 3> center = node->center;
  double width = node->width;
  double offset = 2 * width;
  std::array<double, 3> diff = p - center;

  std::vector<std::array<double, 3>> ret = {center};
  // the cube containing the 8 interpolation nodes is defined by
  // two nodes on the diagonal
  //
  // we already have one of the corner nodes which is the current node
  // find the vector pointing to the other corner
  diff[0] = (diff[0] > 0) ? offset : -offset;
  diff[1] = (diff[1] > 0) ? offset : -offset;
  diff[2] = (diff[2] > 0) ? offset : -offset;

  std::array<double, 3> corner = center + diff;
  ret.push_back(corner);

  // Add the 8 node centers to the centerSet
  for (int i = 0; i < 3; i++) {
    // rotate between axis
    std::array<double, 3> c_copy = std::array<double, 3>(center);
    c_copy[i] += diff[i];
    ret.push_back(c_copy);

    std::array<double, 3> cor_copy = std::array<double, 3>(corner);
    cor_copy[i] -= diff[i];
    ret.push_back(cor_copy);
  }

  return ret;
};

// Compute the neighboring node centers of depth d
std::vector<std::array<double, 3>> nearest_27(Node *node) {
  std::vector<std::array<double, 3>> ret;
  std::array<double, 3> offset = {-2 * node->width, 0, 2 * node->width};

  for (int k = 0; k < 3; k++) {
    for (int j = 0; j < 3; j++) {
      for (int i = 0; i < 3; i++) {
        std::array<double, 3> off = {offset[i], offset[j], offset[k]};
        ret.push_back(node->center + off);
      }
    }
  }

  return ret;
};

// -------------------------------------------------------------------------------------------------//
// Implementation
// -------------------------------------------------------------------------------------------------//

Node::Node(std::vector<id_point> points, std::array<double, 3> center,
           double width, bool is_leaf, int depth)
    : center(center), width(width), depth(depth), is_leaf(is_leaf),
      num_points(points.size()) {
  // children is union type, choose if we want to store points or children
  if (is_leaf) {
    new (&children.points) std::vector<id_point>(points);
  } else {
    new (&children.nodes) std::array<Node *, 8>{nullptr};
  }
};

Node::~Node() {
  if (is_leaf) {
    children.points.~vector();
  } else {
    children.nodes.~array();
  }
};

void Node::Insert(const std::vector<id_point> &points) {
  assert(this->is_leaf);
  this->children.points.insert(this->children.points.end(), points.begin(),
                               points.end());
}

void Node::subdivide() {
  assert(this->is_leaf);
  this->is_leaf = false;

  // get current points
  std::vector<id_point> points = std::move(this->children.points);

  // switch union type
  this->children.points.~vector();
  this->children.nodes = std::array<Node *, 8>{nullptr};

  std::array<std::vector<id_point>, 8> point_partition =
      partition_points(this, points);
  std::vector<std::array<double, 3>> n_centers = split_centers(this);

  // populate children as Node pointers
  for (int i = 0; i < 8; i++) {
    this->children.nodes[i] = new Node(point_partition[i], n_centers[i],
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
  register_node(ret_node);

  if (is_leaf) {
    return ret_node;
  }

  // Figure out which points go in which subdivision
  auto sub_d = partition_points(ret_node, points);

  // Make pointers to subdivision nodes
  auto &n_child = ret_node->children.nodes;
  auto n_centers = split_centers(ret_node);

  for (int i = 0; i < 8; i++) {
    n_child[i] = build(sub_d[i], n_centers[i], width / 2.0, depth + 1);
  };

  return ret_node;
};

Octree::Octree(std::vector<std::array<double, 3>> points, int max_depth,
               int min_depth)
    : _size(points.size()), _points(points), _max_depth(max_depth),
      _min_depth(min_depth) {

  _nodes = std::vector<std::vector<Node *>>(_max_depth + 1);
  _deleted_node_ids = std::vector<std::vector<int>>(_max_depth + 1);

  if (points.size() == 0) {
    this->_root = nullptr;
    return;
  }

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
  double pad = (1.3 * width) / std::pow(2, max_depth - 3);

  this->_root = Octree::build(id_points, center, width / 2 + pad, 0);
};

// priority queue data stores either node or point
// based on distance to query
struct pqData {
  double priority;
  bool is_point;
  int id;

  Node *getNode() { return data.node; };
  std::array<double, 3> getPoint() { return data.pt; };

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

std::vector<int> Octree::kNearestNeighbors(const std::array<double, 3> &query,
                                           int k, const Metric &metric) const {
  std::vector<int> ret;

  std::priority_queue<pqData, std::vector<pqData>, std::greater<pqData>> pq;
  pq.push(pqData(0, this->_root));

  if (k <= 0) {
    k = std::numeric_limits<double>::infinity();
  }

  while (ret.size() < k && !pq.empty()) {

    pqData item = pq.top();
    pq.pop();

    if (item.is_point) {
      // pq item is a point, we add id as a nearest neighbor
      if (!metric(query, item.getPoint())) {
        return ret;
      }
      ret.push_back(item.id);

    } else {
      // add each child of node to pq by distance to query

      // two cases: if node is a leaf or if node is not a leaf
      Node *node = item.getNode();

      if (node->is_leaf) {
        // we add the points of the leaf
        for (const id_point &p : node->children.points) {
          int id = std::get<0>(p);
          std::array<double, 3> cur_p = std::get<1>(p);
          pq.push(pqData(distance(query, cur_p), cur_p, id));
        }
      } else {

        // we add children node to pq
        for (const auto &child : node->children.nodes) {
          if (child != nullptr) {
            pq.push(pqData(distance(query, child), child));
          }
        }
      }
    }
  };

  return ret;
};

std::vector<std::vector<int>>
Octree::kNearestNeighbors(const std::vector<std::array<double, 3>> &queries,
                          int k) const {
  int num_queries = queries.size();
  std::vector<std::vector<int>> res(num_queries);
#pragma omp parallel for
  {
    for (int i = 0; i < num_queries; i++) {
      res[i] = kNearestNeighbors(queries[i], k);
    }
  }

  return res;
};

// Needs to be updated to keep track of _leaf_nodes
int Octree::Delete(Node *node, std::array<double, 3> p) {
  int del_id = -1;
  if (node->is_leaf) {
    std::vector<id_point> &pts = node->children.points;

    for (auto itr = pts.begin(); itr != pts.end(); itr++) {
      id_point id_p = *itr;
      if (p == std::get<1>(id_p)) {
        del_id = std::get<0>(id_p);
        pts.erase(itr);
        node->num_points -= 1;
        break;
      };
    };
    return del_id;
  }

  // find out which child node the point is in
  int idx = node_index_map(node, p);
  Node *next_node = node->children.nodes[idx];
  del_id = Delete(next_node, p);

  if (del_id != -1) {
    node->num_points -= 1;

    if (node->num_points == 0) {
      node->is_leaf = true;
      for (int i = 0; i < 8; i++) {
        // clean up
        delete node->children.nodes[i];
        node->children.nodes[i] = nullptr;
      }

      // change union data type
      node->children.nodes.~array();
    }
  };
  return del_id;
};

void Octree::Delete(std::array<double, 3> p) {
  int del_id = Delete(_root, p);
  if (del_id != -1) {
    _size -= 1;
    _deleted_point_ids.push_back(del_id);
  };
};

void Octree::register_node(Node *node) {
  int depth = node->depth;
  std::vector<Node *> &lvl_nodes = getNodesAtDepth(depth);
  std::vector<int> unused_node_ids = getUnusedNodeIds(depth);

  if (unused_node_ids.size() == 0) {
    node->depth_id = lvl_nodes.size();
    lvl_nodes.push_back(node);
    return;
  }

  node->depth_id = unused_node_ids.back();
  unused_node_ids.pop_back();
  lvl_nodes[node->depth_id] = node;
}

void Octree::unregister_node(Node *node) {
  std::vector<int> unused_node_ids = getUnusedNodeIds(node->depth);
  unused_node_ids.push_back(node->depth_id);
  node->depth_id = -1;
}

int Octree::node_count() const {
  int total;
  for (const auto &node_levels : _nodes) {
    total += node_levels.size();
  }
  return total;
}

struct pqData2 {
  double priority;
  Node *node;

  pqData2(double p, Node *node) : priority(p), node(node) {};

  friend bool operator<(const pqData2 &lhs, const pqData2 &rhs) {
    return lhs.priority < rhs.priority;
  }

  friend bool operator>(const pqData2 &lhs, const pqData2 &rhs) {
    return lhs.priority > rhs.priority;
  }
};

std::vector<int> Octree::RadiusSearch(const std::array<double, 3> &center,
                                      double r, int depth) {
  std::vector<int> found_ids;
  std::priority_queue<pqData2, std::vector<pqData2>, std::greater<pqData2>>
      min_pq;
  min_pq.push(pqData2(0, _root));

  while (!min_pq.empty()) {
    pqData2 data = min_pq.top();
    Node *node = data.node;
    min_pq.pop();

    if (node->is_leaf || node->depth == depth) {
      if (node->depth == depth) {
        found_ids.push_back(node->depth_id);
      }
    } else {
      std::array<Node *, 8> children = node->children.nodes;
      for (int i = 0; i < 8; i++) {
        Node *child = children[i];
        if (child != nullptr) {
          double w = child->width * 1.5;
          std::array<double, 3> c2 = child->center;
          std::array<double, 3> diff = {std::abs(center[0] - c2[0]),
                                        std::abs(center[1] - c2[1]),
                                        std::abs(center[2] - c2[2])};

          double dist =
              std::max(diff[0] - w, std::max(diff[1] - w, diff[2] - w));
          if (dist <= r) {
            min_pq.push(pqData2(dist, child));
          }
        }
      }
    }
  }
  return found_ids;
};

void Octree::Refine(Node *node,
                    const std::vector<std::array<double, 3>> &points,
                    std::vector<Node *> &refinement) {
  if (node->depth_id < 0) {
    this->register_node(node);
  }

  if (points.size() == 0) {
    return;
  }

  if (node->depth == _max_depth) {
    refinement.push_back(node);
    return;
  }

  if (!node->is_leaf || node->depth < _min_depth) {
    bool is_leaf = node->is_leaf;
    if (is_leaf) {
      node->subdivide();
      node->is_leaf = false;
    }

    // Figure out which points go in which subdivision
    auto point_partition = partition_points(node, points);
    std::array<Node *, 8> &children = node->children.nodes;
    std::vector<std::array<double, 3>> n_centers = split_centers(node);

    for (int i = 0; i < 8; i++) {
      if (node->depth == _max_depth - 1 && point_partition[i].size() == 0 &&
          is_leaf) {
        // is child node at max depth
        if (children[i]->children.points.size() == 0) {
          children[i] = nullptr;
        }
        continue;
      }
      if (children[i] == nullptr) {
        children[i] = new Node({}, n_centers[i], node->width / 2.0, true,
                               node->depth + 1);
      }
      Refine(children[i], point_partition[i], refinement);
    }
  }
}

std::vector<Node *>
Octree::Refine(const std::vector<std::array<double, 3>> &points) {
  std::vector<Node *> refinement;

  // Assumes that we do not need to expand Octree bounds
  Refine(_root, points, refinement);
  return refinement;
};

std::vector<Node *> Octree::Neighbors(Node *node) {
  std::vector<std::array<double, 3>> neighbor_c = nearest_27(node);
  std::vector<Node *> neighbors;
  for (int j = 0; j < neighbor_c.size(); j++) {
    Node *found = seek_node(_root, neighbor_c[j], node->depth);
    if (found != nullptr) {
      neighbors.push_back(found);
    }
  }
  return neighbors;
};

std::ostream &operator<<(std::ostream &ofs, const id_point &a) {
  return ofs << std::get<1>(a);
};

std::ostream &operator<<(std::ostream &ofs, const std::array<Node *, 8> &a) {
  if (a[0] == nullptr) {
    ofs << "[null";
  } else {
    ofs << std::endl;
    ofs << *a[0];
  }
  for (int i = 1; i < 8; i++) {
    if (a[i] == nullptr) {
      ofs << ", null";
    } else {
      ofs << std::endl << *a[i];
    }
  };
  if (a[0] == nullptr) {
    ofs << "]";
  }
  return ofs;
};

std::ostream &operator<<(std::ostream &ofs, const Node &n) {
  std::string indent(2 * n.depth, ' ');
  std::string label = (n.is_leaf ? "Leaf" : "Branch");
  ofs << indent << label << "[" << n.num_points << "] {" << std::endl;
  ofs << indent << "  "
      << "center: " << n.center << "," << std::endl;
  ofs << indent << "  "
      << "width: " << n.width << "," << std::endl;
  if (n.is_leaf) {
    ofs << indent << "  " << n.children.points << std::endl;
  } else {
    ofs << indent << "  "
        << "children:" << n.children.nodes << std::endl;
  }
  ofs << indent << "}";
  return ofs;
};

std::ostream &operator<<(std::ostream &ofs, Node *n) {
  std::string indent(2 * n->depth, ' ');
  std::string label = (n->is_leaf ? "Leaf" : "Branch");
  ofs << indent << label << "[" << n->num_points << "] {" << std::endl;
  ofs << indent << "  "
      << "center: " << n->center << "," << std::endl;
  ofs << indent << "  "
      << "width: " << n->width << "," << std::endl;
  if (n->is_leaf) {
    ofs << indent << "  " << n->children.points << std::endl;
  } else {
    ofs << indent << "  "
        << "children:" << n->children.nodes << std::endl;
  }
  ofs << indent << "}";
  return ofs;
};

std::ostream &operator<<(std::ostream &ofs, const Octree &o) {
  ofs << "Octree: [" << std::endl;
  ofs << *o.root() << std::endl;
  ofs << "]";
  return ofs;
}
