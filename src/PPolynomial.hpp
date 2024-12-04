#ifndef PPOLYNOMIAL_HPP
#define PPOLYNOMIAL_HPP

#include "Polynomial.hpp"
#include <iostream>
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

template <int Degree>
bool PPolynomial<Degree>::operator==(const PPolynomial<Degree> &p) {
  return (p.intervals == intervals) && (p.polys == polys);
}

template <int Degree>
bool PPolynomial<Degree>::operator!=(const PPolynomial<Degree> &p) {
  return (p.intervals != intervals) || (p.polys != polys);
}

template <int Degree>
PPolynomial<Degree>
PPolynomial<Degree>::operator+(const PPolynomial<Degree> &p) {
  const std::vector<double> &i1 = intervals;
  const std::vector<double> &i2 = p.intervals;

  const std::vector<Polynomial<Degree>> &p1 = polys;
  const std::vector<Polynomial<Degree>> &p2 = p.polys;

  std::vector<Polynomial<Degree>> new_polys(p1.size() + p2.size());
  std::vector<double> new_intervals;

  int i = i1.size() - 1;
  int j = i2.size() - 1;

  while (i > -1 || j > -1) {
    if (j < 0) {
      new_polys[new_intervals.size()] += p1[i];
      new_intervals.push_back(i1[i--]);
    } else if (i < 0) {
      new_polys[new_intervals.size()] += p2[j];
      new_intervals.push_back(i2[j--]);
    } else {
      new_polys[new_intervals.size()] += p1[i];
      new_polys[new_intervals.size()] += p2[j];
      if (i1[i] > i2[j]) {
        new_intervals.push_back(i1[i--]);
      } else if (i1[i] < i2[j]) {
        new_intervals.push_back(i2[j--]);
      } else {
        new_intervals.push_back(i1[i--]);
        j--;
      }
    }
  }
  std::reverse(new_intervals.begin(), new_intervals.end());
  new_polys.resize(new_intervals.size());
  std::reverse(new_polys.begin(), new_polys.end());

  return PPolynomial<Degree>(new_polys, new_intervals);
}

template <int Degree>
PPolynomial<Degree>
PPolynomial<Degree>::operator-(const PPolynomial<Degree> &p) {
  const std::vector<double> &i1 = intervals;
  const std::vector<double> &i2 = p.intervals;

  const std::vector<Polynomial<Degree>> &p1 = polys;
  const std::vector<Polynomial<Degree>> &p2 = p.polys;

  std::vector<Polynomial<Degree>> new_polys(p1.size() + p2.size());
  std::vector<double> new_intervals;

  int i = i1.size() - 1;
  int j = i2.size() - 1;

  while (i > -1 || j > -1) {
    if (j < 0) {
      new_polys[new_intervals.size()] += p1[i];
      new_intervals.push_back(i1[i--]);
    } else if (i < 0) {
      new_polys[new_intervals.size()] -= p2[j];
      new_intervals.push_back(i2[j--]);
    } else {
      new_polys[new_intervals.size()] += p1[i];
      new_polys[new_intervals.size()] -= p2[j];
      if (i1[i] > i2[j]) {
        new_intervals.push_back(i1[i--]);
      } else if (i1[i] < i2[j]) {
        new_intervals.push_back(i2[j--]);
      } else {
        new_intervals.push_back(i1[i--]);
        j--;
      }
    }
  }
  std::reverse(new_intervals.begin(), new_intervals.end());
  new_polys.resize(new_intervals.size());
  std::reverse(new_polys.begin(), new_polys.end());

  return PPolynomial<Degree>(new_polys, new_intervals);
}

template <int Degree>
std::ostream &operator<<(std::ostream &os, const PPolynomial<Degree> &p) {
  if (p.intervals.size() == 0) {
    return os << "[";
  }
  if (p.intervals.size() == 1) {
    return os << "[ " << p.polys[0] << " { " << p.intervals[0] << " <= x }";
  }

  int end = p.intervals.size() - 1;
  os << "/  " << p.polys[0] << " { " << p.intervals[0]
     << " <= x <= " << p.intervals[1] << " }" << std::endl;
  for (int i = 1; i < end; i++) {
    os << "|  " << p.polys[i] << " { " << p.intervals[i]
       << " <= x <= " << p.intervals[i + 1] << " }" << std::endl;
  }
  return os << "\\  " << p.polys[end] << " { " << p.intervals[end] << " <= x }";
}

#endif
