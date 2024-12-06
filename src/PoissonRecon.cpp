#include "PoissonRecon.hpp"
#include "BSpline.hpp"
#include "HRefine.hpp"
#include "utils/io.hpp"
#include <Eigen/SparseQR>

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
  std::vector<Node *> base_nodes = _octree.getNodesAtDepth(_depth);
  int node_count = base_nodes.size();
  _field_normals = std::vector<std::array<double, 3>>(node_count);
  _field_centers = std::vector<std::array<double, 3>>(node_count);

  for (int i = 0; i < node_count; i++) {

    Node *cur_node = base_nodes[i];
    _field_centers[i] = cur_node->center;

    std::vector<id_point> points = cur_node->children.points;
    for (id_point id_p : points) {
      int p_id = std::get<0>(id_p);
      std::array<double, 3> p = std::get<1>(id_p);

      // get interpolation node centers
      auto n_centers = nearest_8(cur_node, p);

      // get interpolation nodes
      std::vector<Node *> interp_nodes = {cur_node};
      for (int j = 1; j < 8; j++) {
        interp_nodes.push_back(seek_node(_octree.root(), n_centers[j]));
      };

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

        int node_id = node->depth_id;
        _field_normals[node_id] =
            _field_normals[node_id] + weight * _inward_normals[p_id];
      }
    }
  };
}
void PoissonRecon::run() {
  HRefine hr(_octree, _field_normals, BSpline);
  hr.Refine();
  std::vector<double> leaf_coeff = hr.getCoeffAtDepth(_depth);
  _x = Eigen::Map<Eigen::VectorXd>(leaf_coeff.data(), leaf_coeff.size());
  for (int i = 0; i < _depth + 1; i++) {
    std::vector<Node *> cur_nodes = _octree.getNodesAtDepth(i);
    std::vector<std::array<double, 3>> centers;
    for (Node *node : cur_nodes) {
      centers.push_back(node->center);
    }
    _centers.push_back(centers);
  }

  for (int i = 0; i < _depth + 1; i++) {
    save_points(_centers[i],
                "data/centers_depth_" + std::to_string(i) + ".txt");
    writeVectorToFile(hr.coeff[i],
                      "data/x_depth_" + std::to_string(i) + ".txt");
  }
}

//
// void PoissonRecon::run() {
//
//   // get depth d nodes
//   std::vector<Node *> nodes = _octree.getNodesAtDepth(_depth);
//   int node_count = nodes.size();
//
//   _centers = std::vector<std::array<double, 3>>(node_count);
//   _v = Eigen::VectorXd::Zero(node_count);
//
//   std::vector<Eigen::Triplet<double>> triplet_list;
//   triplet_list.reserve(node_count * 8);
//
//   for (int i = 0; i < node_count; i++) {
//     double res = 0;
//     Node *node = nodes[i];
//     std::vector<Node *> neighbors = _octree.Neighbors(node);
//     for (Node *neighbor : neighbors) {
//       res += dot(projection(_divergence_field, node, neighbor),
//                  _field_normals[neighbor->depth_id]);
//       std::array<double, 3> l_p = projection(_laplacian_field, node,
//       neighbor); double entry = l_p[0] + l_p[1] + l_p[2]; if (std::abs(entry)
//       > TOL) {
//         triplet_list.push_back(
//             Eigen::Triplet<double>(node->depth_id, neighbor->depth_id,
//             entry));
//       }
//     };
//     _centers[i] = node->center;
//     _v[i] = res;
//   };
//
//   _L = Eigen::SparseMatrix<double>(node_count, node_count);
//   _L.setFromTriplets(triplet_list.begin(), triplet_list.end());
//
//   // solve the system for x
//   Eigen::ConjugateGradient<Eigen::SparseMatrix<double>,
//                            Eigen::Lower | Eigen::Upper,
//                            Eigen::DiagonalPreconditioner<double>>
//       solver;
//   solver.compute(_L);
//   _x = solver.solve(_v);
// }

void PoissonRecon::write() {
  // int node_count = _centers.size();
  save_points(_points, "points.txt");
  //
  // std::cout << "Node count: " << node_count << std::endl;
  // std::cout << "Matrix size: " << node_count * node_count << std::endl;
  // std::cout << "Num non-zero entries: " << _L.nonZeros() << std::endl;
  //
  // save_sparse_matrix(_L, "L.txt");
  // writeVectorToFile(_v, "v.txt");

  // double max_entry = _x.cwiseAbs().maxCoeff();
  // writeVectorToFile(_x, "x.txt");
  // writeVectorToFile(_x / max_entry, "x_normalized.txt");
}
