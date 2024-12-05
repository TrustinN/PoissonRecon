#ifndef BSPLINE_HPP
#define BSPLINE_HPP

#include "PPolynomial.hpp"

constexpr static double MAX_OUTER_T = 1.5;
constexpr static double MAX_INNER_T = 0.5;
constexpr static double MIN_INNER_T = -0.5;
constexpr static double MIN_OUTER_T = -1.5;

static PPolynomial<2> BSpline = PPolynomial<2>(
    {Polynomial<2>(),
     Polynomial<2>({std::pow(MAX_OUTER_T, 2) / 2.0, MAX_OUTER_T, 0.5}),
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
  std::array<PPolynomial<Degree>, DIM> polys;

  ScalarField() {};
  ScalarField(const ScalarField &sf) : polys(sf.polys) {};
  ScalarField(std::array<PPolynomial<Degree>, DIM> polys) : polys(polys) {};
  ScalarField(const PPolynomial<Degree> &basisF);
  ScalarField(const PPolynomial<Degree> &basisF,
              const std::array<double, DIM> &center);
  ScalarField(const PPolynomial<Degree> &basisF,
              const std::array<double, DIM> &center, double scale);

  ScalarField partialDerivative(int dim) const;
  double integral() const;
  ScalarField<2 * Degree, DIM>
  operator*(const ScalarField<Degree, DIM> &s) const;
  double innerProduct(const ScalarField<Degree, DIM> sf) const;
};

template <int Degree, int DIM>
ScalarField<Degree, DIM>::ScalarField(const PPolynomial<Degree> &bF) {
  for (int i = 0; i < DIM; i++) {
    polys[i] = bF;
  }
};

template <int Degree, int DIM>
ScalarField<Degree, DIM>::ScalarField(const PPolynomial<Degree> &bF,
                                      const std::array<double, DIM> &center) {
  for (int i = 0; i < DIM; i++) {
    polys[i] = bF.shift(center[i]);
  }
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
ScalarField<Degree, DIM>
ScalarField<Degree, DIM>::partialDerivative(int dim) const {
  ScalarField<Degree, DIM> q(*this);
  q.polys[dim] = q.polys[dim].derivative_keep_dim();
  return q;
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
ScalarField<2 * Degree, DIM>
ScalarField<Degree, DIM>::operator*(const ScalarField<Degree, DIM> &sf) const {
  std::array<PPolynomial<2 * Degree>, DIM> new_polys;
  for (int i = 0; i < DIM; i++) {
    new_polys[i] = polys[i] * sf.polys[i];
  }
  return ScalarField<2 * Degree, DIM>(new_polys);
};

template <int Degree, int DIM>
double ScalarField<Degree, DIM>::innerProduct(
    const ScalarField<Degree, DIM> sf) const {
  ScalarField<2 * Degree, DIM> prod = *this * sf;
  return prod.integral();
};

#endif
