#include "pOctree.hpp"
#include "Octree.hpp"
#include "utils/io.hpp"
#include "utils/linalg.hpp"
#include <iostream>

// -------------------------------------------------------------------------------------------------//
// Helper functions
// -------------------------------------------------------------------------------------------------//

Node *seek_node(Node *node, const std::array<double, 3> &p) {

  Node *r_node = node;
  while (!r_node->is_leaf) {
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

pOctree::pOctree(std::vector<std::array<double, 3>> points, int depth)
    : Octree(points, depth, depth) {
  // Also need to refine tree to ensure every point has 8 neighboring nodes
  // at smallest depth

  // defines the center of nodes we need to insert
  std::set<std::array<double, 3>> ins_ctr_set;

  // defines the centers we already have
  std::set<std::array<double, 3>> curr_ctr_set;

  for (Node *node : getNodesAtDepth(_max_depth)) {
    curr_ctr_set.insert(node->center);

    // compute 8 closest center nodes
    for (id_point p : node->children.points) {
      std::vector<std::array<double, 3>> n_centers =
          nearest_8(node, std::get<1>(p));
      ins_ctr_set.insert(n_centers.begin() + 1, n_centers.end());
    };
  }

  // take set difference insert nodes
  std::vector<std::array<double, 3>> diff;
  std::set_difference(ins_ctr_set.begin(), ins_ctr_set.end(),
                      curr_ctr_set.begin(), curr_ctr_set.end(),
                      std::inserter(diff, diff.begin()));

  // insert nodes
  // set refine to be true to have a blank insert
  this->Insert<true>(diff);
};

void pOctree::AssignVecField(std::vector<std::array<double, 3>> normals) {
  // a faster implementation would be to loop over the _base_nodes instead
  for (int i = 0; i < _points.size(); i++) {
    auto p = _points[i];

    // find which node it is in
    Node *cur_node = seek_node(_root, p);

    // get interpolation node centers
    std::vector<std::array<double, 3>> n_centers = nearest_8(cur_node, p);

    // get interpolation nodes
    std::vector<Node *> interp_nodes = {cur_node};
    for (int j = 1; j < 8; j++) {
      interp_nodes.push_back(seek_node(_root, n_centers[j]));
    };

    for (int k = 0; k < 8; k++) {
      int prev_count = _field_centers.size();
      _field_centers.insert(n_centers[k]);
      if (_field_centers.size() > prev_count) {
        _field_nodes.push_back(interp_nodes[k]);
      }
    }

    // trilinear interpolation
    for (Node *node : interp_nodes) {
      // compute distances to center
      std::array<double, 3> diff = p - node->center;
      diff[0] = abs(diff[0]);
      diff[1] = abs(diff[1]);
      diff[2] = abs(diff[2]);

      // invert distance by center distance
      double dist = 2 * node->width;
      double weight = ((dist - diff[0]) / dist) * ((dist - diff[1]) / dist) *
                      ((dist - diff[2]) / dist);

      node->normal = node->normal + weight * normals[i];
    }
  };
};

std::vector<Node *> pOctree::Neighbors(Node *node) {
  std::vector<Node *> ret;
  // start crawling down the tree, when it narrows sufficiently, we can split
  // the search
  Node *cur_node = _root;
  double threshold = 4 * node->width;

  // while (true) {
  //   int idx = node_index_map(cur_node, node->center);
  //   Node *new_node = cur_node->info.children[idx];
  //   std::array<double, 3> diff = new_node->center - node->center;
  //   double dist = std::max(diff[0], std::max(diff[1], diff[2]));
  //   if (new_node->width - dist < threshold) {
  //     break;
  //   }
  //   cur_node = new_node;
  // };

  // start the split search
  std::vector<std::array<double, 3>> targets = nearest_27(node);
  for (int i = 0; i < targets.size(); i++) {
    Node *found = seek_node(cur_node, targets[i]);
    if (found->depth == node->depth) {
      ret.push_back(found);
    }
  }

  return ret;
};
