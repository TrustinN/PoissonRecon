#ifndef POLYNOMIAL_HPP
#define POLYNOMIAL_HPP

#include <array>
#include <cassert>

template <int Degree> struct Polynomial {
  Polynomial() {};
  Polynomial(std::array<double, Degree + 1> coefficients)
      : coefficients(coefficients) {};
  template <int Degree2> Polynomial(Polynomial<Degree2> p);

  std::array<double, Degree + 1> coefficients{0};

  template <int Degree2> bool operator==(const Polynomial<Degree2> &poly) const;
  template <int Degree2> bool operator!=(const Polynomial<Degree2> &poly) const;

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

template <int Degree1>
template <int Degree2>
Polynomial<Degree1>::Polynomial(Polynomial<Degree2> p2) {
  coefficients = std::array<double, Degree1 + 1>{0};
  int max_iter = std::max(Degree1 + 1, Degree2 + 1);
  for (int i = 0; i < max_iter; i++) {
    coefficients[i] = p2.coefficients[i];
  }
}

template <int Degree1>
template <int Degree2>
bool Polynomial<Degree1>::operator==(const Polynomial<Degree2> &poly) const {
  assert(Degree1 == Degree2);
  return poly.coefficients == coefficients;
}

template <int Degree1>
template <int Degree2>
bool Polynomial<Degree1>::operator!=(const Polynomial<Degree2> &poly) const {
  if (Degree1 != Degree2) {
    return true;
  }
  return !(*this == Polynomial<Degree1>(poly));
}

template <int Degree>
Polynomial<Degree>
Polynomial<Degree>::operator+(const Polynomial<Degree> &poly) {
  auto new_coefficients = std::array<double, Degree + 1>(coefficients);
  for (int i = 0; i < Degree + 1; i++) {
    new_coefficients[i] += poly.coefficients[i];
  }

  return Polynomial<Degree>(new_coefficients);
}

template <int Degree>
Polynomial<Degree>
Polynomial<Degree>::operator-(const Polynomial<Degree> &poly) {
  auto new_coefficients = std::array<double, Degree + 1>(coefficients);
  for (int i = 0; i < Degree + 1; i++) {
    new_coefficients[i] -= poly.coefficients[i];
  }

  return Polynomial<Degree>(new_coefficients);
}

template <int Degree>
template <int Degree2>
Polynomial<Degree + Degree2>
Polynomial<Degree>::operator*(const Polynomial<Degree2> &poly) {
  std::array<double, Degree + Degree2 + 1> new_coefficients;
  for (int i = 0; i < Degree + 1; i++) {
    for (int j = 0; j < Degree2 + 1; j++) {
      new_coefficients[i + j] += coefficients[i] * poly.coefficients[j];
    }
  }
  return Polynomial(new_coefficients);
};

template <int Degree>
Polynomial<Degree> &
Polynomial<Degree>::operator+=(const Polynomial<Degree> &poly) {
  for (int i = 0; i < Degree + 1; i++) {
    coefficients[i] += poly.coefficients[i];
  }

  return *this;
}

template <int Degree>
Polynomial<Degree> &
Polynomial<Degree>::operator-=(const Polynomial<Degree> &poly) {
  for (int i = 0; i < Degree + 1; i++) {
    coefficients[i] -= poly.coefficients[i];
  }

  return *this;
}

template <int Degree>
Polynomial<Degree> Polynomial<Degree>::operator+(double s) {
  auto new_coefficients = std::array<double, Degree + 1>(coefficients);
  new_coefficients[0] += s;
  return Polynomial<Degree>(new_coefficients);
}

template <int Degree>
Polynomial<Degree> Polynomial<Degree>::operator-(double s) {
  auto new_coefficients = std::array<double, Degree + 1>(coefficients);
  new_coefficients[0] -= s;
  return Polynomial<Degree>(new_coefficients);
}

template <int Degree>
Polynomial<Degree> Polynomial<Degree>::operator*(double s) {
  auto new_coefficients = std::array<double, Degree + 1>(coefficients);
  for (int i = 0; i < Degree + 1; i++) {
    new_coefficients[i] *= s;
  }
  return Polynomial<Degree>(new_coefficients);
}

template <int Degree>
Polynomial<Degree> Polynomial<Degree>::operator/(double s) {
  auto new_coefficients = std::array<double, Degree + 1>(coefficients);
  for (int i = 0; i < Degree + 1; i++) {
    new_coefficients[i] /= s;
  }
  return Polynomial<Degree>(new_coefficients);
}

template <int Degree>
Polynomial<Degree> &Polynomial<Degree>::operator+=(double s) {
  coefficients[0] += s;
  return *this;
}

template <int Degree>
Polynomial<Degree> &Polynomial<Degree>::operator-=(double s) {
  coefficients[0] -= s;
  return *this;
}

template <int Degree>
Polynomial<Degree> &Polynomial<Degree>::operator*=(double s) {
  for (int i = 0; i < Degree + 1; i++) {
    coefficients[i] *= s;
  }
  return *this;
}

template <int Degree>
Polynomial<Degree> &Polynomial<Degree>::operator/=(double s) {
  for (int i = 0; i < Degree + 1; i++) {
    coefficients[i] /= s;
  }
  return *this;
}

template <int Degree> Polynomial<Degree> Polynomial<Degree>::scale(double s) {
  auto new_coefficients = std::array<double, Degree + 1>(coefficients);
  double f = 1.0;
  for (int i = 0; i < Degree + 1; i++) {
    new_coefficients[i] *= f;
    f /= s;
  }
  return Polynomial<Degree>(new_coefficients);
}

template <int Degree> Polynomial<Degree> Polynomial<Degree>::shift(double t) {
  auto new_coefficients = std::array<double, Degree + 1>(coefficients);
  for (int i = 0; i < Degree + 1; i++) {
    double f = 1.0;
    for (int j = i; j > -1; j--) {
      new_coefficients[j] += f * new_coefficients[i];
      f *= -t * j;
      f /= i - j + 1;
    }
  }
  return Polynomial<Degree>(new_coefficients);
}

template <int Degree> double Polynomial<Degree>::operator()(double x) {
  double res;
  for (int i = Degree; i > -1; i--) {
    res *= x;
    res += coefficients[i];
  }
  return res;
};

template <int Degree> Polynomial<Degree + 1> Polynomial<Degree>::integral() {
  std::array<double, Degree + 1> new_coefficients;
  for (int i = 1; i < Degree + 2; i++) {
    new_coefficients[i] += coefficients[i - 1] / i;
  }
  return Polynomial<Degree + 1>(new_coefficients);
};

template <int Degree> Polynomial<Degree - 1> Polynomial<Degree>::derivative() {
  std::array<double, Degree - 1> new_coefficients;
  for (int i = 0; i < Degree; i++) {
    new_coefficients[i] += i * coefficients[i + 1];
  }
  return Polynomial<Degree - 1>(new_coefficients);
};

template <int Degree> double Polynomial<Degree>::integral(double a, double b) {
  double res = 0.0;
  double c0 = a;
  double c1 = b;
  for (int i = 0; i < Degree + 1; i++) {
    res += coefficients[i] * (c1 - c0) / (i + 1);
    c0 *= a;
    c1 *= b;
  }
  return res;
}

#endif
