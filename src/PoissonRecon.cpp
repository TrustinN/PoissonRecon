#include "PoissonRecon.hpp"
#include <Eigen/SparseQR>

constexpr static double TOL = 1e-6;

PoissonRecon::PoissonRecon(
    const std::vector<std::array<double, 3>> &points,
    const std::vector<std::array<double, 3>> &normals,
    const std::vector<std::array<double, 3>> &inward_normals, int depth)
    : _points(points), _normals(normals), _inward_normals(inward_normals),
      _depth(depth) {

  _octree = pOctree(points, depth);
  _octree.AssignVecField(normals);
};

void PoissonRecon::run() {

  // get depth d nodes
  std::vector<Node *> nodes = _octree.getNodesAtDepth(_depth);
  int node_count = nodes.size();

  _centers = std::vector<std::array<double, 3>>(node_count);
  _v = Eigen::VectorXd::Zero(node_count);

  std::vector<Eigen::Triplet<double>> triplet_list;
  triplet_list.reserve(node_count);
  auto field_normals = _octree.field_normals();

  for (int i = 0; i < node_count; i++) {
    std::cout << i << std::endl;
    double res = 0;
    Node *node = nodes[i];
    std::vector<Node *> neighbors = _octree.Neighbors(node);
    for (Node *neighbor : neighbors) {
      res += dot(projection(_divergence_field, node, neighbor),
                 field_normals[neighbor->depth_id]);
      std::array<double, 3> l_p = projection(_laplacian_field, node, neighbor);
      double entry = l_p[0] + l_p[1] + l_p[2];
      if (std::abs(entry) > TOL) {
        triplet_list.push_back(
            Eigen::Triplet<double>(node->depth_id, neighbor->depth_id, entry));
      }
    };
    _centers[i] = node->center;
    _v[i] = res;
  };

  _L = Eigen::SparseMatrix<double>(node_count, node_count);
  _L.setFromTriplets(triplet_list.begin(), triplet_list.end());

  // solve the system for x
  Eigen::ConjugateGradient<Eigen::SparseMatrix<double>,
                           Eigen::Lower | Eigen::Upper,
                           Eigen::DiagonalPreconditioner<double>>
      solver;
  solver.compute(_L);
  _x = solver.solve(_v);
}

void PoissonRecon::write() {
  int node_count = _centers.size();
  save_points(_centers, "centers.txt");
  save_points(_points, "points.txt");

  std::cout << "Node count: " << node_count << std::endl;
  std::cout << "Matrix size: " << node_count * node_count << std::endl;
  // std::cout << "Num non-zero entries: " << triplet_list.size() << std::endl;

  save_sparse_matrix(_L, "L.txt");
  writeVectorToFile(_v, "v.txt");

  double max_entry = _x.cwiseAbs().maxCoeff();
  writeVectorToFile(_x, "x.txt");
  writeVectorToFile(_x / max_entry, "x_normalized.txt");
}
