#ifndef POLYNOMIALXD_HPP
#define POLYNOMIALXD_HPP

#include "Polynomial.hpp"
#include "utils/linalg.hpp"
#include <numeric>

template <int Base, int Exp> constexpr int int_pow() {
  if constexpr (Exp == 0) {
    return 1;
  } else {
    return Base * int_pow<Base, Exp - 1>();
  }
}

template <int Degree, int DIM = 3> struct PolynomialXD {
  static constexpr int size = int_pow<Degree + 1, DIM>();
  std::array<double, size> coefficients;
  PolynomialXD() { coefficients.fill(0.0); };
  PolynomialXD(const std::array<double, size> &coeff) : coefficients(coeff) {};
  PolynomialXD(const std::array<Polynomial<Degree>, DIM> &p);

  double operator()(const std::array<double, DIM> &p) const;
  PolynomialXD operator+(const PolynomialXD &p) const;
};

template <int Degree, int DIM>
PolynomialXD<Degree, DIM>::PolynomialXD(
    const std::array<Polynomial<Degree>, DIM> &p) {
  for (int d = 0; d < size; d++) {
    int cur = d;
    int ind;
    double coeff = 1.0;
    for (int idx = 0; idx < DIM; idx++) {
      ind = cur % (Degree + 1);
      cur /= (Degree + 1);
      coeff *= p[idx].coefficients[ind];
    }
    coefficients[d] = coeff;
  }
};

// benchmark this maybe
template <int Degree, int DIM>
double
PolynomialXD<Degree, DIM>::operator()(const std::array<double, DIM> &p) const {
  std::array<double, size> res = std::array<double, size>(coefficients);

  std::array<int, DIM> block_size;
  int cur_size = 1.0;
  for (int d = 0; d < DIM; d++) {
    block_size[d] = cur_size;
    cur_size *= (Degree + 1);
  }

  std::array<double, DIM> pows;
  pows.fill(1.0);
  for (int i = 0; i < DIM; i++) {
    int block = block_size[i];
    int stride = block * (Degree + 1);
    for (int d = 0; d < Degree + 1; d++) {
      int start = d * block;
      while (start < size) {
        for (int b = 0; b < block; b++) {
          res[start + b] *= pows[i];
        }
        start += stride;
      }

      pows[i] *= p[i];
    }
  }

  return std::accumulate(res.begin(), res.end(), 0.0);
};

template <int Degree, int DIM>
PolynomialXD<Degree, DIM>
PolynomialXD<Degree, DIM>::operator+(const PolynomialXD<Degree, DIM> &p) const {
  return PolynomialXD<Degree, DIM>(coefficients + p.coefficients);
};

template <int Degree, int DIM>
std::ostream &operator<<(std::ostream &os, const PolynomialXD<Degree, DIM> &p) {
  assert(DIM <= 3);
  char var[3] = {'x', 'y', 'z'};

  if (p.size > 0) {
    os << p.coefficients[0];
  }

  for (int i = 1; i < p.size; i++) {
    int cur = i;
    int ind;
    if (p.coefficients[i] > 0) {
      os << " + " << p.coefficients[i];
    } else if (p.coefficients[i] < 0) {
      os << " - " << std::abs(p.coefficients[i]);
    }
    if (p.coefficients[i] != 0) {
      for (int idx = 0; idx < DIM; idx++) {
        ind = cur % (Degree + 1);
        cur /= (Degree + 1);
        if (ind > 0) {
          os << var[idx];
          if (ind > 1) {
            os << "^" << ind;
          }
        }
      }
    }
  }

  return os;
};

#endif
