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
  PPolynomial(const PPolynomial &p) : polys(p.polys), intervals(p.intervals) {
    assert(intervals.size() == polys.size());
  };

  template <int Degree2>
  PPolynomial(const PPolynomial<Degree2> &p) : intervals(p.intervals) {
    for (const auto &pp : p.polys) {
      polys.push_back(Polynomial<Degree>(pp));
    }
  };

  bool operator==(const PPolynomial &poly) const;
  bool operator!=(const PPolynomial &poly) const;

  PPolynomial operator+(const PPolynomial &poly) const;
  PPolynomial operator-(const PPolynomial &poly) const;
  template <int Degree2>
  PPolynomial<Degree + Degree2>
  operator*(const PPolynomial<Degree2> &poly) const;

  PPolynomial operator+(double s) const;
  PPolynomial operator-(double s) const;
  PPolynomial operator*(double s) const;
  PPolynomial operator/(double s) const;

  PPolynomial &operator+=(double s);
  PPolynomial &operator-=(double s);
  PPolynomial &operator*=(double s);
  PPolynomial &operator/=(double s);

  PPolynomial scale(double s) const;
  PPolynomial shift(double t) const;

  double operator()(double x) const;

  PPolynomial<Degree + 1> integral() const;
  PPolynomial<Degree - 1> derivative() const;
  PPolynomial<Degree> derivative_keep_dim() const;

  double integral(double a, double b) const;
};

template <int Degree>
bool PPolynomial<Degree>::operator==(const PPolynomial<Degree> &p) const {
  return (p.intervals == intervals) && (p.polys == polys);
}

template <int Degree>
bool PPolynomial<Degree>::operator!=(const PPolynomial<Degree> &p) const {
  return (p.intervals != intervals) || (p.polys != polys);
}

template <int Degree>
PPolynomial<Degree>
PPolynomial<Degree>::operator+(const PPolynomial<Degree> &p) const {
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
      new_polys[new_intervals.size()] += p1[i] + p2[j];
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
PPolynomial<Degree>::operator-(const PPolynomial<Degree> &p) const {
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
      new_polys[new_intervals.size()] += p1[i] - p2[j];
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
template <int Degree2>
PPolynomial<Degree + Degree2>
PPolynomial<Degree>::operator*(const PPolynomial<Degree2> &p) const {
  const std::vector<double> &i1 = intervals;
  const std::vector<double> &i2 = p.intervals;

  const std::vector<Polynomial<Degree>> &p1 = polys;
  const std::vector<Polynomial<Degree>> &p2 = p.polys;

  std::vector<Polynomial<Degree + Degree2>> new_polys(p1.size() + p2.size());
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
      new_polys[new_intervals.size()] += p1[i] * p2[j];
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

  return PPolynomial<Degree + Degree2>(new_polys, new_intervals);
};

template <int Degree>
PPolynomial<Degree> PPolynomial<Degree>::operator+(double s) const {
  PPolynomial<Degree> q(*this);
  for (int i = 0; i < polys.size(); i++) {
    q.polys[i] += s;
  }
  return q;
};

template <int Degree>
PPolynomial<Degree> PPolynomial<Degree>::operator-(double s) const {
  PPolynomial<Degree> q(*this);
  for (int i = 0; i < polys.size(); i++) {
    q.polys[i] -= s;
  }
  return q;
};

template <int Degree>
PPolynomial<Degree> PPolynomial<Degree>::operator*(double s) const {
  PPolynomial<Degree> q(*this);
  for (int i = 0; i < polys.size(); i++) {
    q.polys[i] *= s;
  }
  return q;
};

template <int Degree>
PPolynomial<Degree> PPolynomial<Degree>::operator/(double s) const {
  PPolynomial<Degree> q(*this);
  for (int i = 0; i < polys.size(); i++) {
    q.polys[i] /= s;
  }
  return q;
};

template <int Degree>
PPolynomial<Degree> &PPolynomial<Degree>::operator+=(double s) {
  for (int i = 0; i < polys.size(); i++) {
    polys[i] += s;
  }
  return *this;
};

template <int Degree>
PPolynomial<Degree> &PPolynomial<Degree>::operator-=(double s) {
  for (int i = 0; i < polys.size(); i++) {
    polys[i] -= s;
  }
  return *this;
};

template <int Degree>
PPolynomial<Degree> &PPolynomial<Degree>::operator*=(double s) {
  for (int i = 0; i < polys.size(); i++) {
    polys[i] *= s;
  }
  return *this;
};

template <int Degree>
PPolynomial<Degree> &PPolynomial<Degree>::operator/=(double s) {
  for (int i = 0; i < polys.size(); i++) {
    polys[i] /= s;
  }
  return *this;
};

template <int Degree>
PPolynomial<Degree> PPolynomial<Degree>::scale(double s) const {

  PPolynomial<Degree> q(*this);
  for (int i = 0; i < polys.size(); i++) {
    q.polys[i] = polys[i].scale(s);
    q.intervals[i] /= s;
  }
  return q;
};

template <int Degree>
PPolynomial<Degree> PPolynomial<Degree>::shift(double t) const {
  PPolynomial<Degree> q(*this);
  for (int i = 0; i < polys.size(); i++) {
    q.polys[i] = polys[i].shift(t);
    q.intervals[i] += t;
  }
  return q;
};

template <int Degree> double PPolynomial<Degree>::operator()(double x) const {
  for (int i = intervals.size() - 1; i > -1; i--) {
    if (x > intervals[i]) {
      return polys[i](x);
    }
  }
  assert(false);
  return 0;
};

template <int Degree>
PPolynomial<Degree + 1> PPolynomial<Degree>::integral() const {
  std::vector<Polynomial<Degree + 1>> new_polys;
  for (int i = 0; i < polys.size(); i++) {
    new_polys.push_back(polys[i].integral());
  }
  return PPolynomial<Degree + 1>(new_polys, intervals);
};

template <int Degree>
PPolynomial<Degree - 1> PPolynomial<Degree>::derivative() const {
  std::vector<Polynomial<Degree - 1>> new_polys;
  for (int i = 0; i < polys.size(); i++) {
    new_polys.push_back(polys[i].derivative());
  }
  return PPolynomial<Degree - 1>(new_polys, intervals);
};

template <int Degree>
PPolynomial<Degree> PPolynomial<Degree>::derivative_keep_dim() const {
  std::vector<Polynomial<Degree>> new_polys;
  for (int i = 0; i < polys.size(); i++) {
    new_polys.push_back(polys[i].derivative_keep_dim());
  }
  return PPolynomial<Degree>(new_polys, intervals);
};

template <int Degree>
double PPolynomial<Degree>::integral(double a, double b) const {
  int c = 1;
  if (a > b) {
    c = -1;
    double tmp = a;
    a = b;
    b = tmp;
  }
  double res = 0.0;
  for (int i = intervals.size() - 1; i > -1; i--) {

    if (b > intervals[i]) {
      if (a >= intervals[i]) {
        return c * polys[i].integral(a, b);
      }

      while (a < intervals[i]) {
        res += polys[i].integral(intervals[i], b);
        b = intervals[i];
        i--;
      }
      res += polys[i].integral(a, intervals[i + 1]);
      return c * res;
    }
  }

  return c * res;
};

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
