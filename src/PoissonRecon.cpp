#include "PoissonRecon.hpp"
#include "BSpline.hpp"
#include "HRefine.hpp"
#include "Metrics.hpp"
#include "utils/io.hpp"
#include "utils/linalg.hpp"
#include <Eigen/SparseQR>
#include <omp.h>

constexpr static double TOL = 1e-6;

PoissonRecon::PoissonRecon(
    const std::vector<std::array<double, 3>> &points,
    const std::vector<std::array<double, 3>> &normals,
    const std::vector<std::array<double, 3>> &inward_normals, int depth)
    : _points(points), _normals(normals), _inward_normals(inward_normals),
      _depth(depth) {

  _octree = pOctree(points, depth);
  computeVectorField();
};

void PoissonRecon::computeVectorField() {
  std::vector<Node *> &base_nodes = _octree.getNodesAtDepth(_depth);
  int node_count = base_nodes.size();
  _field_normals =
      std::vector<std::array<double, 3>>(node_count, {0.0, 0.0, 0.0});
  std::cout << node_count << std::endl;

  for (int i = 0; i < node_count; i++) {
    Node *cur_node = base_nodes[i];

    std::vector<id_point> points = cur_node->children.points;
    for (id_point id_p : points) {
      int p_id = std::get<0>(id_p);
      std::array<double, 3> p = std::get<1>(id_p);

      // get interpolation node centers
      auto n_centers = nearest_8(cur_node, p);
      int prev_size = base_nodes.size();
      std::vector<Node *> interp_nodes = _octree.Refine(n_centers);
      int new_size = base_nodes.size();
      for (int j = prev_size; j != new_size; j++) {
        _field_normals.push_back({0.0, 0.0, 0.0});
      }

      // trilinear interpolation
      for (Node *node : interp_nodes) {
        // compute distances to center
        std::array<double, 3> diff = p - node->center;

        // invert distance by center distance
        double dist = 2 * node->width;
        double weight = ((dist - std::abs(diff[0])) / dist) *
                        ((dist - std::abs(diff[1])) / dist) *
                        ((dist - std::abs(diff[2])) / dist);

        int node_id = node->depth_id;
        _field_normals[node_id] =
            _field_normals[node_id] + weight * _normals[p_id];
      }
    }
  };

  _field_centers = std::vector<std::array<double, 3>>(base_nodes.size());
  int zero_counter = 0;
  std::array<double, 3> zero = {0.0, 0.0, 0.0};
  for (int i = 0; i < base_nodes.size(); i++) {
    _field_centers[i] = base_nodes[i]->center;
    if (_field_normals[i] == zero) {
      zero_counter++;
    }
  }
  std::cout << _field_centers.size() << std::endl;
  std::cout << zero_counter << std::endl;
  int new_zero_counter = 0;
  for (int i = node_count; i < base_nodes.size(); i++) {
    if (_field_normals[i] == zero) {
      new_zero_counter++;
    }
  }
  std::cout << new_zero_counter << std::endl;
}

void PoissonRecon::run() {
  HRefine hr = HRefine(_octree, _field_normals, BSpline);
  hr.Refine();
  _coeff = hr._coeff;
  _scalar_fields = hr._basis_functions;
  // change back to i < _depth + 1 when we fix addition of PPolynomialXD
  for (int i = 0; i < 4; i++) {
    std::vector<double> cur_coeff = _coeff[i];
    std::vector<ScalarField<2>> cur_fields = _scalar_fields[i];
    std::cout << "Depth: " << i << std::endl;

    for (int j = 0; j < cur_coeff.size(); j++) {
      cur_fields[j] *= cur_coeff[j];
      _indicator_function =
          _indicator_function + PPolynomialXD<2, 3>(cur_fields[j]);
    }
  }
}

void PoissonRecon::write() {
  save_points(_points, "points.txt");

  std::vector<std::vector<double>> widths;
  std::vector<std::vector<double>> iso_vals;
  for (int i = 0; i < _depth + 1; i++) {
    std::vector<Node *> cur_nodes = _octree.getNodesAtDepth(i);

    std::vector<std::array<double, 3>> centers;
    std::vector<double> ws;
    std::vector<double> i_vals;
    for (Node *node : cur_nodes) {
      centers.push_back(node->center);
      ws.push_back(node->width);
      // i_vals.push_back(_indicator_function(node->center));
    }
    _centers.push_back(centers);
    widths.push_back(ws);
    // iso_vals.push_back(i_vals);
  }

  for (int i = 0; i < _depth + 1; i++) {
    save_points(_centers[i],
                "data/centers_depth_" + std::to_string(i) + ".txt");
    writeVectorToFile(_coeff[i], "data/x_depth_" + std::to_string(i) + ".txt");
    writeVectorToFile(widths[i],
                      "data/widths_depth_" + std::to_string(i) + ".txt");
    // writeVectorToFile(iso_vals[i],
    //                   "data/iso_vals_depth_" + std::to_string(i) + ".txt");
  }
}
