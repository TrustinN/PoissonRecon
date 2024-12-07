#include "HRefine.hpp"
#include "utils/io.hpp"
#include "utils/linalg.hpp"
#include <Eigen/Dense>
#include <Eigen/Sparse>

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

  projectRefine(getCoeffAtDepth(depth - 1), getCoeffAtDepth(depth), coarse);
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

  // Compute v
  Eigen::VectorXd v = Eigen::VectorXd::Zero(nodes.size());
  for (int i = 0; i < nodes.size(); i++) {
    Node *cur_node = nodes[i];
    ScalarField<2> cur_node_basisf = cur_node_fields[cur_node->depth_id];
    std::vector<int> v_field_nodes_active =
        _tree.RadiusSearch(cur_node->center, 1.5 * cur_node->width);

    // inner product between gradient of vec field and Node basis
    double v_i = 0.0;
    for (int ii : v_field_nodes_active) {
      Node *v_field_node = v_field_nodes[ii];
      ScalarField<2> v_field_basisf = v_field_fields[v_field_node->depth_id];
      std::array<double, 3> divergence{
          cur_node_basisf.innerProduct(v_field_basisf.partialDerivative(0)),
          cur_node_basisf.innerProduct(v_field_basisf.partialDerivative(1)),
          cur_node_basisf.innerProduct(v_field_basisf.partialDerivative(2))};
      v_i += dot(divergence, _vector_field_normals[ii]);
    }

    v[i] = v_i;
  }

  std::cout << "Computed v!" << std::endl;

  std::vector<Eigen::Triplet<double>> triplet_list;
  triplet_list.reserve(nodes.size() * 8);

  for (int i = 0; i < nodes.size(); i++) {
    Node *cur_node = nodes[i];
    std::vector<Node *> neighbors = _tree.Neighbors(cur_node);

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
  for (int i = 0; i < coarse.size(); i++) {
    Node *parent = coarse[i];
    double contribution = 0.0;
    const std::array<Node *, 8> &children = parent->children.nodes;
    int num_active = 0;

    for (Node *child : children) {
      if (child != nullptr) {
        contribution += fineCoeff[child->depth_id];
        num_active += 1;
      }
    }

    if (num_active > 0) {
      double factor = 1.0 / num_active;
      double residual = coarseCoeff[i] - contribution;
      coarseCoeff[i] -= residual * factor;

      double redistributed_residual = residual * factor;
      for (Node *child : children) {
        if (child != nullptr) {
          fineCoeff[child->depth_id] += redistributed_residual;
        }
      }
    }
  }
}
