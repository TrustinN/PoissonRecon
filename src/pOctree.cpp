#include "pOctree.hpp"
#include "Octree.hpp"
#include "utils/linalg.hpp"
#include <set>

constexpr static double EPSILON = 1e-8;

// -------------------------------------------------------------------------------------------------//
// Helper functions
// -------------------------------------------------------------------------------------------------//

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

Node *pruneNodes(Node *start,
                 const std::vector<std::array<double, 3>> &points) {
  assert(points.size() > 0);
  Node *r_node = start;
  while (true) {
    if (r_node == nullptr) {
      break;
    }
    int idx = node_index_map(r_node, points[0]);
    for (int i = 1; i < points.size(); i++) {
      if (node_index_map(r_node, points[i]) != idx) {
        return r_node;
      };
    }
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

pOctree::pOctree(std::vector<std::array<double, 3>> points, int depth)
    : Octree(points, depth, depth) {
        // Also need to refine tree to ensure every point has 8 neighboring
        // nodes
        // at smallest depth

        // auto max_depth_nodes = getNodesAtDepth(_max_depth);
        // std::vector<std::array<double, 3>> refinement_centers;
        // refinement_centers.reserve(2 * points.size());
        // for (Node *node : max_depth_nodes) {
        //   // compute 8 closest center nodes
        //   for (id_point p : node->children.points) {
        //     std::vector<std::array<double, 3>> n_centers =
        //         nearest_8(node, std::get<1>(p));
        //     for (int i = 0; i < 8; i++) {
        //       refinement_centers.push_back(n_centers[i]);
        //     }
        //   };
        // }
        // this->Insert<true>(refinement_centers);
      };

std::vector<Node *> pOctree::Neighbors(Node *node) {
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
