#ifndef BSPLINE_HPP
#define BSPLINE_HPP

#include "PPolynomial.hpp"

constexpr static double MAX_OUTER_T = 1.5;
constexpr static double MAX_INNER_T = 0.5;
constexpr static double MIN_INNER_T = -0.5;
constexpr static double MIN_OUTER_T = -1.5;

static PPolynomial<2> BSpline = PPolynomial<2>(
    {Polynomial<2>(),
     Polynomial<2>({std::pow(MAX_OUTER_T, 2) / 2.0, -MAX_OUTER_T, 0.5}),
     Polynomial<2>{{MAX_OUTER_T / 2.0, 0, -1.0}},
     Polynomial<2>({std::pow(MAX_OUTER_T, 2) / 2.0, -MAX_OUTER_T, 0.5}),
     Polynomial<2>()},
    {-std::numeric_limits<double>::infinity(), MIN_OUTER_T, MIN_INNER_T,
     MAX_INNER_T, MAX_OUTER_T});

template <int Degree>
inline PPolynomial<Degree> basisFFactory(const PPolynomial<Degree> &p,
                                         double shift, double scale) {
  return p.scale(scale).shift(shift);
};

template <int Degree, int DIM = 3> struct ScalarField {
  std::array<PPolynomial<Degree>, 3> polys{};

  ScalarField() {};
  ScalarField(const PPolynomial<Degree> &basisF,
              const std::array<double, DIM> &center, double scale);

  ScalarField &partialDerivative(int dim);
  double integral() const;
  ScalarField operator*(const ScalarField &s) const;
  double innerProduct(const ScalarField<Degree, DIM> sf) const;
};

template <int Degree, int DIM>
ScalarField<Degree, DIM>::ScalarField(const PPolynomial<Degree> &bF,
                                      const std::array<double, DIM> &center,
                                      double scale) {
  for (int i = 0; i < DIM; i++) {
    polys[i] = basisFFactory(bF, center[i], scale);
  }
};

template <int Degree, int DIM>
ScalarField<Degree, DIM> &ScalarField<Degree, DIM>::partialDerivative(int dim) {
  polys[dim] = polys[dim].derivative_keep_dim();
  return *this;
};

template <int Degree, int DIM>
double ScalarField<Degree, DIM>::integral() const {
  double res = 1.0;
  for (int i = 0; i < DIM; i++) {
    res *= polys[i].integral(-std::numeric_limits<double>::infinity(),
                             std::numeric_limits<double>::infinity());
  }
  return res;
};

template <int Degree, int DIM>
ScalarField<Degree, DIM>
ScalarField<Degree, DIM>::operator*(const ScalarField<Degree, DIM> &sf) const {
  ScalarField<Degree, DIM> res;
  for (int i = 0; i < DIM; i++) {
    res.polys[i] = sf.polys[i] * polys[i];
  }
  return res;
};

template <int Degree, int DIM>
double ScalarField<Degree, DIM>::innerProduct(
    const ScalarField<Degree, DIM> sf) const {
  ScalarField<Degree, DIM> prod = sf * *this;
  return prod.integral();
};

#endif
