#include "HRefine.hpp"
#include "utils/io.hpp"
#include "utils/linalg.hpp"
#include <Eigen/Dense>
#include <Eigen/Sparse>

constexpr static int DEPTH = 2;
constexpr static int MONITOR = 17;
constexpr static bool FULL = false;

double
computeSymmetricError(const Eigen::SparseMatrix<double, Eigen::ColMajor> &L) {
  // Compute the transpose of L
  Eigen::SparseMatrix<double, Eigen::ColMajor> L_transpose = L.transpose();

  // Compute the difference L - L^T
  Eigen::SparseMatrix<double, Eigen::ColMajor> diff = L - L_transpose;

  // Convert to dense to compute Frobenius norm
  Eigen::MatrixXd denseDiff = Eigen::MatrixXd(diff);
  double frobeniusNorm = denseDiff.norm(); // Frobenius norm is the default norm

  return frobeniusNorm;
}

void save_data(Node *cur, Node *leaf, const std::array<double, 3> &normal,
               double weight, const std::string &filename) {
  // store the current Node parameters as functions:
  std::ofstream file(filename, std::ios::app);

  char var[3] = {'x', 'y', 'z'};
  if (file.is_open()) {
    file << "WEIGHT: " << weight << "\n\n";

    for (int i = 0; i < 3; i++) {
      file << "g((x";
      if (cur->center[i] > 0) {
        file << "+";
      };
      file << cur->center[i] << ")/" << cur->width << ")\n";
    }
    file << "\n";

    for (int i = 0; i < 3; i++) {
      file << var[i] << ": " << "g((x";
      if (leaf->center[i] > 0) {
        file << "+";
      };
      file << leaf->center[i] << ")/" << leaf->width << ")\n";
      file << "V_{ector}((" << -leaf->center[i] << ", 0), " << "("
           << 5 * normal[i] - leaf->center[i] << ", 0))";
      file << "\n";
    }
    file << "\n";

    file.close();
    std::cout << "Data saved successfully to " << filename << std::endl;
  } else {
    std::cerr << "Error opening file for writing!" << std::endl;
  }
}

HRefine::HRefine(pOctree tree,
                 const std::vector<std::array<double, 3>> &normals,
                 const PPolynomial<2> &basis)
    : _tree(tree), _basis(basis), _vector_field_normals(normals) {
  _max_depth = tree.max_depth();
  _coeff = std::vector<std::vector<double>>(_max_depth + 1);
  _basis_functions = std::vector<std::vector<ScalarField<2>>>(_max_depth + 1);
  for (int i = 0; i < _max_depth + 1; i++) {
    std::vector<Node *> d_nodes = tree.getNodesAtDepth(i);
    std::vector<ScalarField<2>> &polys = _basis_functions[i];
    double factor = 1.0 / d_nodes[0]->width;
    for (Node *node : d_nodes) {
      polys.push_back(ScalarField<2>(_basis, node->center, factor));
    }
  }
}

void HRefine::Refine(int depth) {
  if (depth > _max_depth) {
    return;
  }

  std::cout << "Refining level " << depth << std::endl;
  const std::vector<Node *> &coarse = _tree.getNodesAtDepth(depth - 1);
  const std::vector<Node *> &fine = _tree.getNodesAtDepth(depth);
  std::vector<double> &coarseCoeff = getCoeffAtDepth(depth - 1);
  std::vector<double> fineCoeff =
      coarseToFineRefine(coarseCoeff, coarse, fine, depth);
  setCoeffAtDepth(fineCoeff, depth);

  Refine(depth + 1);

  projectRefine(getCoeffAtDepth(depth - 1), getCoeffAtDepth(depth), coarse,
                depth);
}

void HRefine::Refine() {
  int depths = _tree.max_depth() + 1;
  _coeff[0] = {1};
  setCoeffAtDepth(computeCoeff(getCoeffAtDepth(0), _tree.getNodesAtDepth(0), 0),
                  0);
  Refine(1);
};

std::vector<double>
HRefine::coarseToFineRefine(std::vector<double> &coarseCoeff,
                            const std::vector<Node *> &coarse,
                            const std::vector<Node *> &fine, int depth) {
  std::vector<double> fineCoeff = initializeRefine(coarseCoeff, coarse, fine);
  fineCoeff = computeCoeff(fineCoeff, fine, depth);
  return fineCoeff;
};

std::vector<double>
HRefine::initializeRefine(const std::vector<double> &coarseCoeff,
                          const std::vector<Node *> &coarse,
                          const std::vector<Node *> &fine) {

  std::vector<double> fineCoeff(fine.size());
  for (int i = 0; i < coarse.size(); i++) {
    Node *parent = coarse[i];
    double cc = coarseCoeff[i] / 8.0;
    for (Node *child : parent->children.nodes) {
      if (child != nullptr) {
        fineCoeff[child->depth_id] = cc;
      }
    }
  }

  return fineCoeff;
};

std::vector<double> HRefine::computeCoeff(std::vector<double> &start,
                                          const std::vector<Node *> &nodes,
                                          int depth) {

  std::vector<Node *> v_field_nodes = _tree.getNodesAtDepth(_max_depth);
  std::vector<ScalarField<2>> v_field_fields = getFieldsAtDepth(_max_depth);
  std::vector<ScalarField<2>> cur_node_fields = getFieldsAtDepth(depth);

  std::vector<std::array<double, 3>> debug_points;
  std::vector<std::array<double, 3>> debug_normals;
  std::vector<std::array<double, 3>> debug_centers;
  std::vector<double> debug_weights;

  // Compute v
  Eigen::VectorXd v = Eigen::VectorXd::Zero(nodes.size());
  for (int i = 0; i < nodes.size(); i++) {
    Node *cur_node = nodes[i];
    if (depth == DEPTH && (i == MONITOR || FULL)) {
      debug_centers.push_back(cur_node->center);
    }
    ScalarField<2> cur_node_basisf = cur_node_fields[cur_node->depth_id];
    std::vector<int> v_field_nodes_active =
        _tree.RadiusSearch(cur_node->center, 1.5 * cur_node->width, _max_depth);

    // inner product between gradient of vec field and Node basis
    double v_i = 0.0;
    for (int ii : v_field_nodes_active) {
      Node *v_field_node = v_field_nodes[ii];
      assert(v_field_node->depth_id == ii);
      ScalarField<2> v_field_basisf = v_field_fields[v_field_node->depth_id];
      std::array<double, 3> divergence{
          cur_node_basisf.innerProduct(v_field_basisf.partialDerivative(0)),
          cur_node_basisf.innerProduct(v_field_basisf.partialDerivative(1)),
          cur_node_basisf.innerProduct(v_field_basisf.partialDerivative(2))};

      double val = dot(divergence, _vector_field_normals[ii]);
      v_i += val;
      if (depth == DEPTH) {
        if (i == MONITOR || FULL) {
          debug_points.push_back(v_field_node->center);
          debug_normals.push_back(_vector_field_normals[ii]);
          // debug_weights.push_back(abs(divergence[0]) + abs(divergence[1]) +
          //                         abs(divergence[2]));
          debug_weights.push_back(val);
        }
      }
    }

    v[i] = v_i;
    if (depth == DEPTH) {
      if (i == MONITOR) {
        std::cout << "v[" << MONITOR << "]: " << v_i << std::endl;
      }
    }
  }

  std::cout << "Computed v!" << std::endl;

  std::vector<Eigen::Triplet<double>> triplet_list;
  triplet_list.reserve(nodes.size() * 8);

  std::vector<int> nonzeros;
  std::vector<double> nonzerosLs;
  for (int i = 0; i < nodes.size(); i++) {
    Node *cur_node = nodes[i];
    // this is causing the weird coefficient calculations
    std::vector<Node *> neighbors = _tree.Neighbors(cur_node);
    // std::vector<Node *> neighbors = {cur_node};

    // inner product of node and neighboring nodes
    ScalarField<2> cur_node_basisf = cur_node_fields[cur_node->depth_id];
    ScalarField<2> cnb_dx =
        cur_node_basisf.partialDerivative(0).partialDerivative(0);
    ScalarField<2> cnb_dy =
        cur_node_basisf.partialDerivative(1).partialDerivative(1);
    ScalarField<2> cnb_dz =
        cur_node_basisf.partialDerivative(2).partialDerivative(2);

    for (Node *neighbor : neighbors) {
      ScalarField<2> neighbor_basisf = cur_node_fields[neighbor->depth_id];
      double L_ij = cnb_dx.innerProduct(neighbor_basisf) +
                    cnb_dy.innerProduct(neighbor_basisf) +
                    cnb_dz.innerProduct(neighbor_basisf);
      if (depth == DEPTH && i == MONITOR) {
        std::cout << "L[" << cur_node->depth_id << "]" << "["
                  << neighbor->depth_id << "]" << " = " << L_ij << std::endl;
        nonzeros.push_back(neighbor->depth_id);
        nonzerosLs.push_back(L_ij);
      }
      triplet_list.push_back(
          Eigen::Triplet<double>(cur_node->depth_id, neighbor->depth_id, L_ij));
    }
  }

  // Compute L
  Eigen::SparseMatrix<double, Eigen::ColMajor> L(nodes.size(), nodes.size());
  L.setFromTriplets(triplet_list.begin(), triplet_list.end());
  std::cout << "Computed L!" << std::endl;
  double error = computeSymmetricError(L);
  std::cout << "Error: " << error << std::endl;
  if (error > 1e-4) {
    std::cout << "Depth of Error: " << depth << std::endl;
    throw std::runtime_error("Symmetry error detected");
  }

  // Solve for x
  Eigen::ConjugateGradient<Eigen::SparseMatrix<double>,
                           Eigen::Lower | Eigen::Upper>
      solver;
  solver.compute(L);
  Eigen::VectorXd guess =
      Eigen::Map<Eigen::VectorXd>(start.data(), start.size());
  Eigen::VectorXd res = solver.solveWithGuess(v, guess);
  std::cout << "Solved x!" << std::endl;

  if (depth == DEPTH) {
    for (int k = 0; k < nodes.size(); k++) {
      Node *cur_node = nodes[k];
      if (res(k) > 0 &&
          distance(cur_node->center, {0, 0, 0}) > std::pow(20, 2)) {
        std::cout << "Monitor: " << k << std::endl;
      }
    }
    std::cout << "x[" << MONITOR << "] = " << res(MONITOR) << std::endl;
    double actual = 0.0;
    for (int p = 0; p < nonzeros.size(); p++) {
      std::cout << "x[" << nonzeros[p] << "] = " << res(nonzeros[p])
                << std::endl;
      actual += res(nonzeros[p]) * nonzerosLs[p];
    }
    std::cout << "v_actual[" << MONITOR << "] = " << actual << std::endl;

    save_points(debug_points, "debug/points.txt");
    save_points(debug_normals, "debug/normals.txt");
    save_points(debug_centers, "debug/centers.txt");
    writeVectorToFile(debug_weights, "debug/weights.txt");
    if (FULL) {
      writeVectorToFile(res, "debug/coeff.txt");
    } else {
      writeVectorToFile({res(MONITOR)}, "debug/coeff.txt");
    }
  }

  return std::vector<double>(res.data(), res.data() + res.size());
};

void HRefine::projectRefine(std::vector<double> &coarseCoeff,
                            std::vector<double> &fineCoeff,
                            const std::vector<Node *> &coarse, int depth) {

  // double factor = 1.0 / 8.0;
  // std::vector<ScalarField<2>> coarse_fields = getFieldsAtDepth(depth - 1);
  // std::vector<ScalarField<2>> fine_fields = getFieldsAtDepth(depth);
  // for (int i = 0; i < coarse.size(); i++) {
  //   Node *parent = coarse[i];
  //   double contribution = 0.0;
  //   const std::array<Node *, 8> &children = parent->children.nodes;
  //
  //   ScalarField<2> parent_field = coarse_fields[parent->depth_id];
  //   ScalarField<2> pdx =
  //   parent_field.partialDerivative(0).partialDerivative(0); ScalarField<2>
  //   pdy = parent_field.partialDerivative(1).partialDerivative(1);
  //   ScalarField<2> pdz =
  //   parent_field.partialDerivative(2).partialDerivative(2);
  //
  //   for (Node *child : children) {
  //     if (child != nullptr) {
  //       ScalarField<2> child_field = fine_fields[child->depth_id];
  //       double inner = pdx.innerProduct(child_field) +
  //                      pdy.innerProduct(child_field) +
  //                      pdz.innerProduct(child_field);
  //       contribution += fineCoeff[child->depth_id] * inner;
  //       fineCoeff[child->depth_id] -=
  //           coarseCoeff[parent->depth_id] * inner * factor;
  //     }
  //   }
  //   coarseCoeff[i] -= contribution;
  // }
}

// void HRefine::projectRefine(std::vector<double> &coarseCoeff,
//                             std::vector<double> &fineCoeff,
//                             const std::vector<Node *> &coarse, int depth) {
//
//   double factor = 1.0 / 8.0;
//   std::vector<ScalarField<2>> coarse_fields = getFieldsAtDepth(depth - 1);
//   std::vector<ScalarField<2>> fine_fields = getFieldsAtDepth(depth);
//   for (int i = 0; i < coarse.size(); i++) {
//     Node *parent = coarse[i];
//     double contribution = 0.0;
//     const std::array<Node *, 8> &children = parent->children.nodes;
//
//     ScalarField<2> parent_field = coarse_fields[parent->depth_id];
//     ScalarField<2> pdx =
//     parent_field.partialDerivative(0).partialDerivative(0); ScalarField<2>
//     pdy = parent_field.partialDerivative(1).partialDerivative(1);
//     ScalarField<2> pdz =
//     parent_field.partialDerivative(2).partialDerivative(2); for (Node *child
//     : children) {
//       if (child != nullptr) {
//         ScalarField<2> child_field = fine_fields[child->depth_id];
//         double inner = pdx.innerProduct(child_field) +
//                        pdy.innerProduct(child_field) +
//                        pdz.innerProduct(child_field);
//         std::cout << inner << std::endl;
//         contribution += fineCoeff[child->depth_id] * inner;
//       }
//     }
//     coarseCoeff[i] -= contribution * factor;
//
//     double redistributed_residual = contribution * factor * factor;
//     for (Node *child : children) {
//       if (child != nullptr) {
//         fineCoeff[child->depth_id] += redistributed_residual;
//       }
//     }
//   }
// }
