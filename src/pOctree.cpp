#include "pOctree.hpp"
#include "utils/linalg.hpp"

#include <set>

pOctree::pOctree(std::vector<std::array<double, 3>> points, int depth)
    : Octree(points, depth, depth) {
  // Also need to refine tree to ensure every point has 8 neighboring nodes
  // at smallest depth

  // defines the center of nodes we need to insert
  std::set<std::array<double, 3>> ins_ctr_set;

  // defines the centers we already have
  std::set<std::array<double, 3>> curr_ctr_set;

  for (Node *node : _base_nodes) {
    std::array<double, 3> center = node->center;
    curr_ctr_set.insert(center);
    double width = node->width;
    double offset = 2 * width;

    // compute 8 closest center nodes
    for (id_point p : node->info.points) {

      // the cube containing the 8 interpolation nodes is defined by
      // two nodes on the diagonal
      //
      // we already have one of the corner nodes which is the current node
      // find the vector pointing to the other corner

      std::array<double, 3> diff = std::get<1>(p) - center;
      diff[0] = (diff[0] > 0) ? offset : -offset;
      diff[1] = (diff[1] > 0) ? offset : -offset;
      diff[2] = (diff[2] > 0) ? offset : -offset;

      std::array<double, 3> corner = center + diff;
      ins_ctr_set.insert(corner);

      // Add the 8 node centers to the centerSet
      for (int i = 0; i < 3; i++) {
        // rotate between axis
        std::array<double, 3> c_copy = std::array<double, 3>(center);
        c_copy[i] += diff[i];
        ins_ctr_set.insert(c_copy);

        std::array<double, 3> cor_copy = std::array<double, 3>(corner);
        cor_copy[i] -= diff[i];
        ins_ctr_set.insert(cor_copy);
      }
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
