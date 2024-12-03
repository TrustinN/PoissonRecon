#include "PPolynomial.hpp"

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

  const std::vector<double> &p1 = polys;
  const std::vector<double> &p2 = p.polys;

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

  const std::vector<double> &p1 = polys;
  const std::vector<double> &p2 = p.polys;

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
