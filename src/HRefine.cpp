#include "HRefine.hpp"
#include "Eigen/Sparse"
#include "utils/linalg.hpp"

HRefine::HRefine(pOctree tree,
                 const std::vector<std::array<double, 3>> &normals,
                 const PPolynomial<2> &basis)
    : tree(tree), basis(basis), _vector_field_normals(normals) {
  coeff = std::vector<std::vector<double>>(tree.max_depth() + 1);
}

void HRefine::Refine() {
  int depths = tree.max_depth() + 1;
  for (int i = 1; i < depths; i++) {
    std::vector<Node *> coarser_nodes = tree.getNodesAtDepth(i - 1);
    std::vector<Node *> cur_nodes = tree.getNodesAtDepth(i);
    coarseToFineRefine(getCoeffAtDepth(i - 1), coarser_nodes, cur_nodes);
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

  std::vector<double> fineCoeff(fine.size(), 0.0);
  for (int i = 0; i < coarse.size(); i++) {
    Node *parent = coarse[i];
    double cc = coarseCoeff[i] / 8.0;
    for (Node *children : parent->children.nodes) {
      fineCoeff[children->depth_id] = cc;
    }
  }

  return fineCoeff;
};

std::vector<double> HRefine::computeCoeff(const std::vector<double> &start,
                                          const std::vector<Node *> &nodes) {

  int max_depth = tree.max_depth();
  std::vector<Node *> v_field_nodes = tree.getNodesAtDepth(max_depth);

  // Compute v
  Eigen::VectorXd v = Eigen::VectorXd::Zero(nodes.size());
  for (int i = 0; i < nodes.size(); i++) {
    Node *cur_node = nodes[i];
    ScalarField<2> cur_node_basisf(basis, cur_node->center, cur_node->width);
    std::vector<int> v_field_nodes_active =
        tree.RadiusSearch(cur_node->center, 1.5 * cur_node->width);

    // inner product between gradient of vec field and Node basis
    double v_i = 0.0;
    for (int ii : v_field_nodes_active) {
      Node *v_field_node = v_field_nodes[ii];
      ScalarField<2> v_field_basisf(basis, v_field_node->center,
                                    v_field_node->width);
      std::array<double, 3> divergence;
      for (int i = 0; i < 3; i++) {
        ScalarField<2> dv_basisf = v_field_basisf.partialDerivative(i);
        divergence[i] = dv_basisf.innerProduct(cur_node_basisf);
      }
      v_i += dot(divergence, _vector_field_normals[ii]);
    }

    v[i] = v_i;
  }

  std::vector<Eigen::Triplet<double>> triplet_list;
  triplet_list.reserve(nodes.size() * 8);

  for (int i = 0; i < nodes.size(); i++) {
    Node *cur_node = nodes[i];
    std::vector<std::array<double, 3>> neighbor_c = nearest_27(cur_node);
    std::vector<Node *> neighbors;
    for (int i = 0; i < neighbor_c.size(); i++) {
      Node *found = seek_node(cur_node, neighbor_c[i], cur_node->depth);
      neighbors.push_back(found);
    }

    // inner product of node and neighboring nodes
    ScalarField<2> cur_node_basisf(basis, cur_node->center, cur_node->width);
    for (Node *neighbor : neighbors) {
      ScalarField<2> neighbor_basisf(basis, neighbor->center, neighbor->width);
      double L_ij = 0.0;
      for (int i = 0; i < 3; i++) {
        L_ij += cur_node_basisf.innerProduct(
            neighbor_basisf.partialDerivative(i).partialDerivative(i));
      }
      triplet_list.push_back(
          Eigen::Triplet<double>(cur_node->depth_id, neighbor->depth_id, L_ij));
    }
  }

  // Compute L
  Eigen::SparseMatrix<double, Eigen::ColMajor> L(nodes.size(), nodes.size());
  L.setFromTriplets(triplet_list.begin(), triplet_list.end());

  // Solve for x
  Eigen::ConjugateGradient<Eigen::SparseMatrix<double>,
                           Eigen::Lower | Eigen::Upper>
      solver;
  solver.compute(L);
  Eigen::VectorXd res = solver.solveWithGuess(v, start);

  return std::vector<double>(res.data(), res.data() + res.size());
};

void HRefine::projectRefine(std::vector<double> &coarseCoeff,
                            std::vector<double> &fineCoeff,
                            const std::vector<Node *> &coarse) {
  double factor = 1.0 / 8.0;
  for (int i = 0; i < coarse.size(); i++) {
    Node *parent = coarse[i];
    double contribution = 0.0;
    const std::array<Node *, 8> &children = parent->children.nodes;
    for (Node *child : children) {
      contribution += fineCoeff[child->depth_id];
    }
    coarseCoeff[i] -= contribution * factor;
    double res = coarseCoeff[i] * factor;
    for (Node *child : children) {
      fineCoeff[child->depth_id] -= res;
    }
  }
};
