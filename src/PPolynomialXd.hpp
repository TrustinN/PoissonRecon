#ifndef PPOLYNOMIALXD_HPP
#define PPOLYNOMIALXD_HPP

#include "BSpline.hpp"
#include "PPolynomial.hpp"
#include "PolynomialXd.hpp"
#include "utils/io.hpp"

template <int Degree, int DIM> struct PPolynomialXD {
  std::array<std::vector<double>, DIM> _intervals;
  std::vector<PolynomialXD<Degree, DIM>> _polys;
  std::array<int, DIM> _sizes;
  std::array<double, DIM + 1> _bit_map_pow{1.0};
  int _divisions = 1;

  PPolynomialXD();
  PPolynomialXD(std::vector<PolynomialXD<Degree, DIM>> polys,
                std::array<std::vector<double>, DIM> intervals);
  PPolynomialXD(const std::array<PPolynomial<Degree>, DIM> &ppolys);
  PPolynomialXD(const ScalarField<2, DIM> &sf) : PPolynomialXD(sf.polys) {};

  double operator()(const std::array<double, DIM> &p) const;

  int getPolyFromIntervalIndices(const std::array<int, DIM> &ind) const;

  PPolynomialXD operator+(const PPolynomialXD &p) const;
};

template <int Degree, int DIM> PPolynomialXD<Degree, DIM>::PPolynomialXD() {
  _polys.push_back(PolynomialXD<Degree, DIM>());
  _sizes.fill(1);
  _intervals.fill({-std::numeric_limits<double>::infinity()});
};

template <int Degree, int DIM>
PPolynomialXD<Degree, DIM>::PPolynomialXD(
    std::vector<PolynomialXD<Degree, DIM>> polys,
    std::array<std::vector<double>, DIM> intervals)
    : _polys(polys), _intervals(intervals) {
  for (int i = 0; i < DIM; i++) {
    _sizes[i] = intervals[i].size();
    _bit_map_pow[i + 1] = _bit_map_pow[i] * _sizes[i];
    _divisions *= _sizes[i];
  }
};

template <int Degree, int DIM>
PPolynomialXD<Degree, DIM>::PPolynomialXD(
    const std::array<PPolynomial<Degree>, DIM> &p) {
  for (int i = 0; i < DIM; i++) {
    _intervals[i] = p[i].intervals;
    _sizes[i] = _intervals[i].size();
    _bit_map_pow[i + 1] = _bit_map_pow[i] * _sizes[i];
    _divisions *= _sizes[i];
  }

  std::array<int, DIM> idxmap;
  idxmap.fill(0);
  for (int d = 0; d < _divisions; d++) {
    std::array<Polynomial<Degree>, DIM> new_polys;

    for (int i = 0; i < DIM; i++) {
      new_polys[i] = p[i].polys[idxmap[i]];
    }
    _polys.push_back(PolynomialXD<Degree, DIM>(new_polys));

    int cur_carry = 0;
    idxmap[cur_carry] += 1;
    while (idxmap[cur_carry] == _sizes[cur_carry]) {
      idxmap[cur_carry] = 0;
      cur_carry += 1;
      idxmap[cur_carry] += 1;
    }
  }
}

template <int Degree, int DIM>
PPolynomialXD<Degree, DIM> PPolynomialXD<Degree, DIM>::operator+(
    const PPolynomialXD<Degree, DIM> &p) const {
  std::array<std::vector<double>, DIM> new_intervals;
  std::array<std::vector<std::pair<int, int>>, DIM> ij;

  // Construct the new intervals first
  for (int d = 0; d < DIM; d++) {
    const std::vector<double> &i1 = _intervals[d];
    const std::vector<double> &i2 = p._intervals[d];
    std::vector<double> &ni = new_intervals[d];
    std::vector<std::pair<int, int>> &pairings = ij[d];

    int i = 0, j = 0;
    while (i < i1.size() || j < i2.size()) {
      if (i >= i1.size()) {
        pairings.push_back({i - 1, j});
        ni.push_back(i2[j++]);
      } else if (j >= i2.size()) {
        pairings.push_back({i, j - 1});
        ni.push_back(i1[i++]);
      } else {
        if (i1[i] < i2[j]) {
          ni.push_back(i1[i++]);
        } else if (i1[i] > i2[j]) {
          ni.push_back(i2[j++]);
        } else {
          ni.push_back(i1[i++]);
          j++;
        }
        pairings.push_back({i - 1, j - 1});
      }
    }
  }

  int divisions = 1;
  std::array<int, DIM> sizes;
  for (int i = 0; i < new_intervals.size(); i++) {
    const auto &vec = new_intervals[i];
    divisions *= vec.size();
    sizes[i] = vec.size();
  }
  std::cout << divisions << std::endl;

  // construct the polys based on the new intervals
  std::vector<PolynomialXD<Degree, DIM>> new_polys;
  std::array<int, DIM> idxmap;
  idxmap.fill(0);

  for (int d = 0; d < divisions; d++) {

    // make two indexing arrays to get the polynomials
    // from the lhs PPolynomialXD and another for the
    // rhs PPolynomialXD
    std::array<int, DIM> indexing_arr1;
    std::array<int, DIM> indexing_arr2;
    for (int i = 0; i < DIM; i++) {
      int idx = idxmap[i];
      // contains the pair (i, j) where i indexes into rhs interval and
      // j indexes into lhs interval
      std::pair<int, int> pair = ij[i][idx];
      indexing_arr1[i] = std::get<0>(pair);
      indexing_arr2[i] = std::get<1>(pair);
    }

    // recover the polynomials from both piecewise inputs
    PolynomialXD<Degree, DIM> p1 =
        _polys[getPolyFromIntervalIndices(indexing_arr1)];
    PolynomialXD<Degree, DIM> p2 =
        p._polys[p.getPolyFromIntervalIndices(indexing_arr2)];
    PolynomialXD<Degree, DIM> p_new = p1 + p2;
    new_polys.push_back(p_new);

    int cur_carry = 0;
    idxmap[cur_carry] += 1;
    while (idxmap[cur_carry] == sizes[cur_carry]) {
      idxmap[cur_carry] = 0;
      cur_carry += 1;
      idxmap[cur_carry] += 1;
    }
  }

  return PPolynomialXD<Degree, DIM>(new_polys, new_intervals);
};

// Remember to replace all evaluations with a binary search
template <int Degree, int DIM>
double
PPolynomialXD<Degree, DIM>::operator()(const std::array<double, DIM> &p) const {
  std::array<int, DIM> ind;
  for (int d = 0; d < DIM; d++) {
    std::vector<double> intv = _intervals[d];
    for (int i = intv.size() - 1; i > -1; i--) {
      if (p[d] > intv[i]) {
        ind[d] = i;
        break;
      }
    }
  }
  int idx = getPolyFromIntervalIndices(ind);
  return _polys[idx](p);
};

template <int Degree, int DIM>
int PPolynomialXD<Degree, DIM>::getPolyFromIntervalIndices(
    const std::array<int, DIM> &ind) const {
  int idx = 0;

  for (int i = 0; i < DIM; i++) {
    idx += ind[i] * _bit_map_pow[i];
  }
  return idx;
};

template <int Degree, int DIM>
std::ostream &operator<<(std::ostream &os,
                         const PPolynomialXD<Degree, DIM> &p) {
  assert(DIM <= 3);
  char vars[3] = {'x', 'y', 'z'};

  std::array<int, DIM> idxmap;
  idxmap.fill(0);
  for (int d = 0; d < p._divisions; d++) {
    os << std::endl << p._polys[d];
    for (int i = 0; i < DIM; i++) {
      int cur_i = idxmap[i];
      std::vector<double> interval = p._intervals[i];
      if (cur_i == p._sizes[i] - 1) {
        os << " { " << interval[cur_i] << " <= " << vars[i] << " }";
      } else {
        os << " { " << interval[cur_i] << " <= " << vars[i]
           << " <= " << interval[cur_i + 1] << " }";
      }
    }

    int cur_carry = 0;
    idxmap[cur_carry] += 1;
    while (idxmap[cur_carry] == p._sizes[cur_carry]) {
      idxmap[cur_carry] = 0;
      cur_carry += 1;
      idxmap[cur_carry] += 1;
    }
  }
  return os;
}

#endif
