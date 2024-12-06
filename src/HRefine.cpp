#include "HRefine.hpp"
#include "BSpline.hpp"
#include "utils/io.hpp"
#include "utils/linalg.hpp"
#include <Eigen/Dense>
#include <Eigen/Sparse>

HRefine::HRefine(pOctree tree,
                 const std::vector<std::array<double, 3>> &normals,
                 const PPolynomial<2> &basis)
    : tree(tree), basis(basis), _vector_field_normals(normals) {
  coeff = std::vector<std::vector<double>>(tree.max_depth() + 1);
}

void HRefine::Refine() {
  int depths = tree.max_depth() + 1;
  coeff[0] = {-1};
  setCoeffAtDepth(computeCoeff(getCoeffAtDepth(0), tree.getNodesAtDepth(0)), 0);
  for (int i = 1; i < depths; i++) {
    std::vector<Node *> coarser_nodes = tree.getNodesAtDepth(i - 1);
    std::vector<Node *> cur_nodes = tree.getNodesAtDepth(i);
    std::cout << "Refining level " << i << std::endl;
    setCoeffAtDepth(
        coarseToFineRefine(getCoeffAtDepth(i - 1), coarser_nodes, cur_nodes),
        i);
  }
};

std::vector<double>
HRefine::coarseToFineRefine(std::vector<double> &coarseCoeff,
                            const std::vector<Node *> &coarse,
                            const std::vector<Node *> &fine) {
  std::vector<double> fineCoeff = initializeRefine(coarseCoeff, coarse, fine);
  fineCoeff = computeCoeff(fineCoeff, fine);
  projectRefine(coarseCoeff, fineCoeff, coarse);
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
                                          const std::vector<Node *> &nodes) {

  int max_depth = tree.max_depth();
  std::vector<Node *> v_field_nodes = tree.getNodesAtDepth(max_depth);

  // Compute v
  Eigen::VectorXd v = Eigen::VectorXd::Zero(nodes.size());
  for (int i = 0; i < nodes.size(); i++) {
    Node *cur_node = nodes[i];
    ScalarField<2> cur_node_basisf(basis, cur_node->center,
                                   1 / cur_node->width);
    std::vector<int> v_field_nodes_active =
        tree.RadiusSearch(cur_node->center, 1.5 * cur_node->width);

    // inner product between gradient of vec field and Node basis
    double v_i = 0.0;
    for (int ii : v_field_nodes_active) {
      Node *v_field_node = v_field_nodes[ii];
      ScalarField<2> v_field_basisf(basis, v_field_node->center,
                                    1 / v_field_node->width);
      std::array<double, 3> divergence{
          cur_node_basisf.innerProduct(v_field_basisf.partialDerivative(0)),
          cur_node_basisf.innerProduct(v_field_basisf.partialDerivative(1)),
          cur_node_basisf.innerProduct(v_field_basisf.partialDerivative(2))};
      v_i += dot(divergence, _vector_field_normals[ii]);
    }

    v[i] = 20000 * v_i;
  }

  std::cout << "Computed v!" << std::endl;

  std::vector<Eigen::Triplet<double>> triplet_list;
  triplet_list.reserve(nodes.size() * 8);

  for (int i = 0; i < nodes.size(); i++) {
    Node *cur_node = nodes[i];
    std::vector<std::array<double, 3>> neighbor_c = nearest_27(cur_node);
    std::vector<Node *> neighbors;
    for (int j = 0; j < neighbor_c.size(); j++) {
      Node *found = seek_node(cur_node, neighbor_c[j], cur_node->depth);
      if (found != nullptr) {
        neighbors.push_back(found);
      }
    }

    // inner product of node and neighboring nodes
    ScalarField<2> cur_node_basisf(basis, cur_node->center,
                                   1 / cur_node->width);
    ScalarField<2> cnb_dx =
        cur_node_basisf.partialDerivative(0).partialDerivative(0);
    ScalarField<2> cnb_dy =
        cur_node_basisf.partialDerivative(1).partialDerivative(1);
    ScalarField<2> cnb_dz =
        cur_node_basisf.partialDerivative(2).partialDerivative(2);

    for (Node *neighbor : neighbors) {
      ScalarField<2> neighbor_basisf(basis, neighbor->center,
                                     1 / neighbor->width);
      double L_ij = cnb_dx.innerProduct(neighbor_basisf) +
                    cnb_dy.innerProduct(neighbor_basisf) +
                    cnb_dz.innerProduct(neighbor_basisf);
      triplet_list.push_back(
          Eigen::Triplet<double>(cur_node->depth_id, neighbor->depth_id, L_ij));
    }
  }

  // Compute L
  Eigen::SparseMatrix<double, Eigen::ColMajor> L(nodes.size(), nodes.size());
  L.setFromTriplets(triplet_list.begin(), triplet_list.end());
  std::cout << "Computed L!" << std::endl;

  // Solve for x
  Eigen::ConjugateGradient<Eigen::SparseMatrix<double>,
                           Eigen::Lower | Eigen::Upper>
      solver;
  solver.compute(L);
  Eigen::VectorXd guess =
      Eigen::Map<Eigen::VectorXd>(start.data(), start.size());
  Eigen::VectorXd res = solver.solveWithGuess(v, guess);
  std::cout << "Solved x!" << std::endl;

  return std::vector<double>(res.data(), res.data() + res.size());
};

void HRefine::projectRefine(std::vector<double> &coarseCoeff,
                            std::vector<double> &fineCoeff,
                            const std::vector<Node *> &coarse) {
  // double factor = 1.0 / 8.0;
  // for (int i = 0; i < coarse.size(); i++) {
  //   Node *parent = coarse[i];
  //   double contribution = 0.0;
  //   const std::array<Node *, 8> &children = parent->children.nodes;
  //   for (Node *child : children) {
  //     if (child != nullptr) {
  //       contribution += fineCoeff[child->depth_id];
  //     }
  //   }
  //   std::cout << "Error: " << std::abs(contribution - coarseCoeff[i])
  //             << std::endl;
  //   coarseCoeff[i] -= contribution * factor;
  //   double res = coarseCoeff[i] * factor;
  //   for (Node *child : children) {
  //     if (child != nullptr) {
  //       fineCoeff[child->depth_id] -= res;
  //     }
  //   }
  // }
};
