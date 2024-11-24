#include "PoissonRecon.hpp"
#include <Eigen/SparseQR>
#include <fstream>
#include <iostream>

void writeVectorToFile(const Eigen::VectorXd &vec,
                       const std::string &filename) {
  std::ofstream file(filename);
  if (!file.is_open()) {
    std::cerr << "Failed to open file: " << filename << std::endl;
    return;
  }

  if (vec.size() > 0) {
    file << vec[0];
    for (int i = 1; i < vec.size(); ++i) {
      file << ", " << vec[i];
    }
  }

  file.close();
  if (file.good()) {
    std::cout << "Vector successfully written to " << filename << std::endl;
  } else {
    std::cerr << "Error occurred while writing to the file." << std::endl;
  }
}
PoissonRecon::PoissonRecon(
    const std::vector<std::array<double, 3>> &points,
    const std::vector<std::array<double, 3>> &normals,
    const std::vector<std::array<double, 3>> &inward_normals, int depth)
    : _points(points), _normals(normals), _inward_normals(inward_normals),
      _depth(depth) {

  _octree = pOctree(points, depth);
  _octree.AssignVecField(normals);

  std::vector<Node *> nodes = _octree.leaf_nodes();
  int node_count = nodes.size();

  _v = Eigen::VectorXd::Zero(node_count);

  std::vector<Eigen::Triplet<double>> triplet_list;
  triplet_list.reserve(node_count);

  for (int i = 0; i < node_count; i++) {
    double res = 0;
    Node *node = nodes[i];
    std::vector<Node *> neighbors = _octree.Neighbors(node);
    for (Node *neighbor : neighbors) {
      res += projection(_divergence_field, node, neighbor);
      double entry = projection(_laplacian_field, node, neighbor) / node->width;
      triplet_list.push_back(
          Eigen::Triplet<double>(node->id, neighbor->id, entry));
    };
    _v[i] = res;
  };
  writeVectorToFile(_v, "v.txt");

  _L = Eigen::SparseMatrix<double>(node_count, node_count);
  _L.setFromTriplets(triplet_list.begin(), triplet_list.end());

  // solve the system for x
  Eigen::ConjugateGradient<Eigen::SparseMatrix<double>,
                           Eigen::Lower | Eigen::Upper>
      solver;
  solver.compute(_L);

  _x = solver.solve(_v);
  writeVectorToFile(_x, "x.txt");
};
