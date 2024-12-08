#ifndef HREFINE_HPP
#define HREFINE_HPP

#include "BSpline.hpp"
#include "PPolynomial.hpp"
#include "pOctree.hpp"
#include <vector>

struct HRefine {
  pOctree _tree;
  int _max_depth;
  PPolynomial<2> _basis;
  std::vector<std::vector<ScalarField<2>>> _basis_functions;
  std::vector<std::array<double, 3>> _vector_field_normals;
  std::vector<std::vector<double>> _coeff;

  HRefine(pOctree tree, const std::vector<std::array<double, 3>> &normals,
          const PPolynomial<2> &basis);

  void setCoeffAtDepth(std::vector<double> dCoeff, int depth) {
    _coeff[depth] = dCoeff;
  };
  std::vector<double> &getCoeffAtDepth(int depth) { return _coeff[depth]; };
  const std::vector<ScalarField<2>> &getFieldsAtDepth(int depth) const {
    return _basis_functions[depth];
  };
  std::vector<double> computeCoeff(std::vector<double> &start,
                                   const std::vector<Node *> &nodes, int depth);

  std::vector<double> initializeRefine(const std::vector<double> &coarseCoeff,
                                       const std::vector<Node *> &coarse,
                                       const std::vector<Node *> &fine);

  std::vector<double> coarseToFineRefine(std::vector<double> &coarseCoeff,
                                         const std::vector<Node *> &coarse,
                                         const std::vector<Node *> &fine,
                                         int depth);

  void projectRefine(std::vector<double> &coarseCoeff,
                     std::vector<double> &fineCoeff,
                     const std::vector<Node *> &coarse, int depth);

  void Refine();
  void Refine(int depth);
};

#endif
