#ifndef BASIS_HPP
#define BASIS_HPP

#include "integration.hpp"
#include <iostream>

struct basisF1 {
  static double operator()(double t);
};

struct dbasisF1 {
  static double operator()(double t);
};

struct d2basisF1 {
  static double operator()(double t);
};

struct divergenceField {
  std::array<double, 27> x_field =
      inner_product_field<basisF1, dbasisF1, basisF1, basisF1, basisF1,
                          basisF1>();
  std::array<double, 27> y_field =
      inner_product_field<basisF1, basisF1, basisF1, dbasisF1, basisF1,
                          basisF1>();
  std::array<double, 27> z_field =
      inner_product_field<basisF1, basisF1, basisF1, basisF1, basisF1,
                          dbasisF1>();
};

struct laplacianField {
  std::array<double, 27> x_field =
      inner_product_field<basisF1, d2basisF1, basisF1, basisF1, basisF1,
                          basisF1>();
  std::array<double, 27> y_field =
      inner_product_field<basisF1, basisF1, basisF1, d2basisF1, basisF1,
                          basisF1>();
  std::array<double, 27> z_field =
      inner_product_field<basisF1, basisF1, basisF1, basisF1, basisF1,
                          d2basisF1>();
};

template <typename field_type>
std::array<double, 3> projection(field_type field, Node *n1, Node *n2) {
  std::array<double, 3> center = n2->center;

  std::array<double, 3> diff = center - n1->center;

  std::array<double, 3> bit_map{0};
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

  double infl_x = field.x_field[idx];
  double infl_y = field.y_field[idx];
  double infl_z = field.z_field[idx];

  return {infl_x, infl_y, infl_z};
};

template <typename field_type>
std::array<double, 3> projection(field_type field, std::array<double, 3> c1,
                                 std::array<double, 3> c2) {
  // std::array<double, 3> center = n2->center;

  std::array<double, 3> diff = c2 - c1;

  std::array<double, 3> bit_map{0};
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

  double infl_x = field.x_field[idx];
  double infl_y = field.y_field[idx];
  double infl_z = field.z_field[idx];

  return {infl_x, infl_y, infl_z};
};

#endif
