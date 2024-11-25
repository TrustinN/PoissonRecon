#include "PoissonRecon.hpp"
#include <Eigen/SparseQR>
#include <fstream>
#include <iostream>

constexpr static double TOL = 1e-6;

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

void save_sparse_matrix(const Eigen::SparseMatrix<double> &matrix,
                        const std::string &filename) {
  // Using the Matrix Market format to save the sparse matrix
  std::ofstream file(filename);
  if (file.is_open()) {
    Eigen::SparseMatrix<double> matrix_transposed =
        matrix.transpose(); // Transpose to make saving easier
    matrix_transposed
        .makeCompressed();     // Ensure the matrix is in compressed form
    file << matrix_transposed; // Write the matrix to file
    file.close();
    std::cout << "Sparse matrix saved successfully to " << filename
              << std::endl;
  } else {
    std::cerr << "Error opening file for writing!" << std::endl;
  }
}

void save_points(const std::vector<std::array<double, 3>> &points,
                 const std::string &filename) {
  std::ofstream file(filename);

  if (file.is_open()) {
    for (const auto &point : points) {
      file << point[0] << " " << point[1] << " " << point[2] << "\n";
    }
    file.close();
    std::cout << "Points saved successfully to " << filename << std::endl;
  } else {
    std::cerr << "Error opening file for writing!" << std::endl;
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

  // get depth d nodes
  std::vector<Node *> leaf_nodes = _octree.leaf_nodes();
  std::vector<Node *> nodes;
  for (int i = 0; i < leaf_nodes.size(); i++) {
    if (leaf_nodes[i]->depth == depth) {
      leaf_nodes[i]->depth_id = nodes.size();
      nodes.push_back(leaf_nodes[i]);
    }
  }

  int node_count = nodes.size();

  std::vector<std::array<double, 3>> centers(node_count);
  _v = Eigen::VectorXd::Zero(node_count);

  std::vector<Eigen::Triplet<double>> triplet_list;
  triplet_list.reserve(node_count);

  for (int i = 0; i < node_count; i++) {
    double res = 0;
    Node *node = nodes[i];
    std::vector<Node *> neighbors = _octree.Neighbors(node);
    for (Node *neighbor : neighbors) {
      res +=
          dot(projection(_divergence_field, node, neighbor), neighbor->normal);
      std::array<double, 3> l_p = projection(_laplacian_field, node, neighbor);
      double entry = (l_p[0] + l_p[1] + l_p[2]) / node->width;
      if (std::abs(entry) > TOL) {
        triplet_list.push_back(
            Eigen::Triplet<double>(node->depth_id, neighbor->depth_id, entry));
      }
    };
    centers[i] = node->center;
    _v[i] = res;
  };
  writeVectorToFile(_v, "v.txt");

  _L = Eigen::SparseMatrix<double>(node_count, node_count);
  _L.setFromTriplets(triplet_list.begin(), triplet_list.end());

  Eigen::MatrixXd denseL = Eigen::MatrixXd(_L);
  Eigen::MatrixXd diff = denseL - denseL.transpose();

  std::cout << "Symmetry error: " << diff.cwiseAbs().maxCoeff() << std::endl;
  std::cout << "Node count: " << node_count << std::endl;
  std::cout << "Matrix size: " << node_count * node_count << std::endl;
  std::cout << "Num non-zero entries: " << triplet_list.size() << std::endl;

  save_sparse_matrix(_L, "L.txt");

  // solve the system for x
  Eigen::ConjugateGradient<Eigen::SparseMatrix<double>,
                           Eigen::Lower | Eigen::Upper,
                           Eigen::DiagonalPreconditioner<double>>
      solver;
  solver.compute(_L);

  _x = solver.solve(_v);
  writeVectorToFile(_x, "x.txt");

  double max_entry = _x.cwiseAbs().minCoeff();
  writeVectorToFile(_x / max_entry, "x_normalized.txt");

  save_points(centers, "centers.txt");
  save_points(points, "points.txt");
};
