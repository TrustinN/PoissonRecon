#ifndef PPOLYNOMIAL_HPP
#define PPOLYNOMIAL_HPP

#include "Polynomial.hpp"
#include <vector>

template <int Degree> struct PPolynomial {
  std::vector<Polynomial<Degree>> polys;
  std::vector<double> intervals;

  PPolynomial() {};
  PPolynomial(std::vector<Polynomial<Degree>> polys,
              std::vector<double> intervals)
      : polys(polys), intervals(intervals) {
    assert(intervals.size() == polys.size());
  };

  bool operator==(const PPolynomial &poly);
  bool operator!=(const PPolynomial &poly);

  PPolynomial operator+(const PPolynomial &poly);
  PPolynomial operator-(const PPolynomial &poly);
  template <int Degree2>
  PPolynomial<Degree + Degree2> operator*(const PPolynomial<Degree2> &poly);

  PPolynomial &operator+=(const PPolynomial &poly);
  PPolynomial &operator-=(const PPolynomial &poly);

  PPolynomial operator+(double s);
  PPolynomial operator-(double s);
  PPolynomial operator*(double s);
  PPolynomial operator/(double s);

  PPolynomial &operator+=(double s);
  PPolynomial &operator-=(double s);
  PPolynomial &operator*=(double s);
  PPolynomial &operator/=(double s);

  PPolynomial scale(double s);
  PPolynomial shift(double t);

  double operator()(double x);

  PPolynomial<Degree + 1> integral();
  PPolynomial<Degree - 1> derivative();

  double integral(double a, double b);
};

#endif
