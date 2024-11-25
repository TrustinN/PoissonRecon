#ifndef BASIS_HPP
#define BASIS_HPP

#include "Octree.hpp"
#include "utils/linalg.hpp"
#include <array>

// 1-D basis function
struct basisF1 {

  // gives the integral over t of dF1(t)/dt * F1(t)
  static double div_weight(double t);
  // gives the integral over t of F1^2(t)
  // computed as the complement to the current integral
  static double div_weight_cmpl(double t);

  // gives d^2F1(t)/dt^2 * F1(t) (Laplacian)
  static double laplace_weight(double t);

  // change of variables from x -> t given node o with center oc and width ow
  static double cov(double x, double oc, double ow) { return (x - oc) / ow; };
  static std::array<double, 3> cov(const std::array<double, 3> &p, double oc,
                                   double ow) {
    return {cov(p[0], oc, ow), cov(p[1], oc, ow), cov(p[2], oc, ow)};
  };

  // Evaluation at p
  static double operator()(double p);
};

// F(x, y, z) = F1(x) * F1(y) * F1(z)
struct basisF {
  double operator()(const std::array<double, 3> &p);
};

struct Field {
  static const std::array<std::vector<int>, 27> loc_to_v_field_idx;
  std::array<double, 3> dw{0};
  std::array<double, 3> wc{0};

  std::array<double, 27> int_field_x{0};
  std::array<double, 27> int_field_y{0};
  std::array<double, 27> int_field_z{0};

  std::array<double, 27> _infl_x{0};
  std::array<double, 27> _infl_y{0};
  std::array<double, 27> _infl_z{0};
};

template <typename field_type> void initialize_field(field_type &f) {
  for (int k = 0; k < 3; k++) {
    for (int j = 0; j < 3; j++) {
      for (int i = 0; i < 3; i++) {
        int idx = i + 3 * j + 9 * k;
        f.int_field_x[idx] = f.dw[i] * f.wc[j] * f.wc[k];
        f.int_field_y[idx] = f.dw[j] * f.wc[k] * f.wc[i];
        f.int_field_z[idx] = f.dw[k] * f.wc[i] * f.wc[j];
      }
    }
  }

  for (int i = 0; i < 27; i++) {
    std::vector<int> infl_areas = f.loc_to_v_field_idx[i];
    for (int pos : infl_areas) {
      f._infl_x[i] += f.int_field_x[pos];
      f._infl_y[i] += f.int_field_y[pos];
      f._infl_z[i] += f.int_field_z[pos];
    }
  }
}

// We output the integral of dF1(t_x)/dt_x * F1(t_x) * (F1(t_y) * F1(t_z))^2 wrt
// t_x, t_y, t_z
//
// array is divided into sections of 3 x 3 x 3 = 27
// on intervals [-1.5, -0.5], [-0.5, 0.5], [0.5, 1.5]
struct divVField : public Field {
  divVField();
};

struct laplaceField : public Field {
  laplaceField();
};

template <typename field_type>
std::array<double, 3> projection(field_type field, Node *n1, Node *n2) {
  std::array<double, 3> center = n2->center;

  std::array<double, 3> diff = center - n1->center;

  std::array<double, 3> bit_map;
  if (diff[0] == 0) {
    bit_map[0] = 1;
  } else if (diff[0] > 0) {
    bit_map[0] = 2;
  }
  if (diff[1] == 0) {
    bit_map[1] = 1;
  } else if (diff[1] > 0) {
    bit_map[1] = 2;
  }
  if (diff[2] == 0) {
    bit_map[2] = 1;
  } else if (diff[2] > 0) {
    bit_map[2] = 2;
  }

  int idx = bit_map[0] + 3 * bit_map[1] + 9 * bit_map[2];

  double infl_x = field._infl_x[idx];
  double infl_y = field._infl_y[idx];
  double infl_z = field._infl_z[idx];

  return {infl_x, infl_y, infl_z};
};

#endif
