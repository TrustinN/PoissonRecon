#ifndef BASIS_HPP
#define BASIS_HPP

#include "Octree.hpp"
#include "utils/linalg.hpp"
#include <array>

struct divergenceField {
  std::array<double, 27> x_field;
  std::array<double, 27> y_field;
  std::array<double, 27> z_field;
  divergenceField();
};

struct laplacianField {
  std::array<double, 27> x_field;
  std::array<double, 27> y_field;
  std::array<double, 27> z_field;
  laplacianField();
};

double integral_f1_f2(double t);
double integral_df1_f2(double t);
double integral_d2f1_f2(double t);

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

#endif
