#ifndef POLYNOMIAL_HPP
#define POLYNOMIAL_HPP

#include <array>

template <int Degree> struct Polynomial {
  Polynomial() {};
  Polynomial(std::array<double, Degree + 1> coefficients)
      : coefficients(coefficients) {};

  std::array<double, Degree + 1> coefficients;

  bool operator==(const Polynomial &poly);
  bool operator!=(const Polynomial &poly);

  Polynomial operator+(const Polynomial &poly);
  Polynomial operator-(const Polynomial &poly);
  template <int Degree2>
  Polynomial<Degree + Degree2> operator*(const Polynomial<Degree2> &poly);

  Polynomial &operator+=(const Polynomial &poly);
  Polynomial &operator-=(const Polynomial &poly);

  Polynomial operator+(double s);
  Polynomial operator-(double s);
  Polynomial operator*(double s);
  Polynomial operator/(double s);

  Polynomial &operator+=(double s);
  Polynomial &operator-=(double s);
  Polynomial &operator*=(double s);
  Polynomial &operator/=(double s);

  Polynomial scale(double s);
  Polynomial shift(double t);

  double operator()(double x);

  Polynomial<Degree + 1> integral();
  Polynomial<Degree - 1> derivative();

  double integral(double a, double b);
};

#endif
