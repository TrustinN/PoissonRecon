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
  int getPolyFromIntervalIndices(const std::array<int, DIM> &ind,
                                 const std::array<double, DIM + 1> &map) const;

  PPolynomialXD operator+(const PPolynomialXD &p) const;
  PPolynomialXD &operator+=(const PPolynomialXD &p);
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

template <typename T> int binary_insert(const std::vector<T> &v, T d) {
  int a = 0, b = v.size() - 1;
  if (d > v[b]) {
    return b;
  }
  while (true) {
    int mid = (a + b) / 2;
    if (d < v[mid]) {
      if (d >= v[mid - 1]) {
        return mid - 1;
      }
      b = mid;
    } else {
      a = mid + 1;
    }
  }
}

// template <int Degree, int DIM>
// PPolynomialXD<Degree, DIM> &
// PPolynomialXD<Degree, DIM>::operator+=(const PPolynomialXD<Degree, DIM> &p) {
//
//   std::array<std::vector<std::pair<int, int>>, DIM> ij;
//
//   // Modify intervals
//   for (int d = 0; d < DIM; d++) {
//     std::vector<double> &i1 = _intervals[d];
//     const std::vector<double> &i2 = p._intervals[d];
//     std::vector<std::pair<int, int>> &pairings = ij[d];
//
//     // This only works for BSplines, we ignore the 0 parts
//     for (int i = 1; i < i2.size() - 1; i++) {
//       int insert_idx = binary_insert<double>(i1, i2[i]);
//       i1.insert(i1.begin() + insert_idx, i2[i]);
//       pairings.push_back({insert_idx, i});
//     }
//   }
//
//   for (auto c : ij) {
//     std::cout << c << std::endl;
//   }
//
//   int divisions = 1;
//   std::array<int, DIM> sizes;
//   std::array<int, DIM> ij_sizes;
//   std::array<double, DIM + 1> bit_map_pow{1.0};
//   for (int i = 0; i < DIM; i++) {
//     int s = _intervals[i].size();
//     divisions *= s;
//     sizes[i] = s;
//     bit_map_pow[i + 1] = bit_map_pow[i] * sizes[i];
//     ij_sizes[i] = ij[i].size();
//   }
//
//   // Modify the poly_list
//   std::vector<PolynomialXD<Degree, DIM>> new_polys;
//   std::vector<int> insert_indexes;
//   // cannot insert immediately as getPolyFromIntervalIndices relies
//   // on unchanged _divisions, _sizes, so on
//   std::array<int, DIM> idxmap;
//   idxmap.fill(0);
//   bool finish = false;
//   while (!finish) {
//     // make two indexing arrays to get the polynomials
//     // from the lhs PPolynomialXD and another for the
//     // rhs PPolynomialXD
//     std::array<int, DIM> indexing_arr1;
//     std::array<int, DIM> indexing_arr2;
//     std::array<int, DIM> ins_index;
//     for (int i = 0; i < DIM; i++) {
//       int cur_idx_of_ij = idxmap[i];
//       std::pair<int, int> pairings = ij[i][cur_idx_of_ij];
//       // Since we inserted into our interval in ascending order
//       ins_index[i] = std::get<0>(pairings);
//
//       indexing_arr1[i] = std::get<0>(pairings) - 1;
//       indexing_arr2[i] = std::get<1>(pairings);
//     }
//
//     // recover the polynomials from both piecewise inputs
//     PolynomialXD<Degree, DIM> p1 =
//         _polys[getPolyFromIntervalIndices(indexing_arr1)];
//     PolynomialXD<Degree, DIM> p2 =
//         p._polys[p.getPolyFromIntervalIndices(indexing_arr2)];
//     PolynomialXD<Degree, DIM> p_new = p1 + p2;
//     new_polys.push_back(p_new);
//     insert_indexes.push_back(
//         getPolyFromIntervalIndices(ins_index, bit_map_pow));
//
//     int cur_carry = 0;
//     idxmap[cur_carry] += 1;
//     while (idxmap[cur_carry] == ij_sizes[cur_carry]) {
//       idxmap[cur_carry] = 0;
//       cur_carry += 1;
//       if (cur_carry == DIM) {
//         finish = true;
//         break;
//       }
//       idxmap[cur_carry] += 1;
//     }
//   }
//
//   std::cout << "Inserting polys" << std::endl;
//   std::cout << "size: " << _polys.size() << std::endl;
//   std::cout << insert_indexes << std::endl;
//   for (int i = 0; i < insert_indexes.size(); i++) {
//     int idx = insert_indexes[i] + i;
//     std::cout << idx << std::endl;
//     _polys.insert(_polys.begin() + idx, new_polys[i]);
//   }
//
//   // Update all other variables
//   _sizes = sizes;
//   _divisions = divisions;
//   _bit_map_pow = bit_map_pow;
//
//   return *this;
// };

// rather than explicit interval / poly management which would be mem
// inefficient store the other interval and polys and perform lazy evaluation
// when needed

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
    ind[d] = binary_insert(intv, p[d]);
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
int PPolynomialXD<Degree, DIM>::getPolyFromIntervalIndices(
    const std::array<int, DIM> &ind,
    const std::array<double, DIM + 1> &map) const {

  int idx = 0;
  for (int i = 0; i < DIM; i++) {
    idx += ind[i] * map[i];
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
