#ifndef HREFINE_HPP
#define HREFINE_HPP

#include "PPolynomial.hpp"
#include "pOctree.hpp"
#include <vector>

struct HRefine {
  pOctree tree;
  PPolynomial<2> basis;
  std::vector<std::array<double, 3>> _vector_field_normals;
  std::vector<std::vector<double>> coeff;

  HRefine(pOctree tree, const std::vector<std::array<double, 3>> &normals,
          const PPolynomial<2> &basis);

  void setCoeffAtDepth(std::vector<double> dCoeff, int depth) {
    coeff[depth] = dCoeff;
  };
  std::vector<double> &getCoeffAtDepth(int depth) { return coeff[depth]; };

  std::vector<double> computeCoeff(std::vector<double> &start,
                                   const std::vector<Node *> &nodes);

  std::vector<double> initializeRefine(const std::vector<double> &coarseCoeff,
                                       const std::vector<Node *> &coarse,
                                       const std::vector<Node *> &fine);

  std::vector<double> coarseToFineRefine(std::vector<double> &coarseCoeff,
                                         const std::vector<Node *> &coarse,
                                         const std::vector<Node *> &fine);

  void projectRefine(std::vector<double> &coarseCoeff,
                     std::vector<double> &fineCoeff,
                     const std::vector<Node *> &coarse);

  void Refine();
};

#endif
