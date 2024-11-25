#include "basis.hpp"
#include <cassert>

constexpr static double MAX_OUTER_T = 1.5;
constexpr static double MAX_INNER_T = 0.5;
constexpr static double MIN_INNER_T = -0.5;
constexpr static double MIN_OUTER_T = -1.5;
constexpr static double EPSILON = 1e-10;

double sign(double x) { return (x > double(0)) - (x < double(0)); };

double basisF1::div_weight(double t) {
  double abs_t = abs(t);

  if (abs_t > MAX_OUTER_T) {
    return 0;

  } else if (abs_t > MAX_INNER_T) {
    return std::pow(MAX_OUTER_T - abs_t, 4) / 8;

  } else {
    return (-MAX_OUTER_T * std::pow(t, 2) + std::pow(t, 4)) / 2;
  };
};

double basisF1::div_weight_cmpl(double t) {
  double abs_t = std::abs(t);

  if (abs_t > MAX_OUTER_T) {
    return 0;

  } else if (abs_t > MAX_INNER_T) {
    return -sign(t) * std::pow(MAX_OUTER_T - abs_t, 5) / 20;

  } else {
    return 0.5625 * t - std::pow(t, 3) / 2 + std::pow(t, 5) / 5;
  };
};

double basisF1::laplace_weight(double t) {
  double abs_t = std::abs(t);
  if (abs_t > MAX_OUTER_T) {
    return 0;

  } else if (abs_t > MAX_INNER_T) {
    return -sign(t) * std::pow(MAX_OUTER_T - abs_t, 3) / 6;

  } else {
    return 2 * std::pow(t, 3) / 3 - MAX_OUTER_T * t;
  }
}

double basisF1::operator()(double t) {
  double abs_t = std::abs(t);

  if (abs_t > MAX_OUTER_T) {
    return 0;

  } else if (abs_t > MAX_INNER_T) {
    return std::pow(MAX_OUTER_T - abs_t, 2) / 2;

  } else {
    return MAX_OUTER_T / 2 - std::pow(t, 2);
  };
}

double basisF::operator()(const std::array<double, 3> &p) {
  return basisF1::operator()(p[0]) * basisF1::operator()(p[1]) *
         basisF1::operator()(p[2]);
};

const std::array<std::vector<int>, 27> Field::loc_to_v_field_idx = {
    // first level
    std::vector<int>{0},
    std::vector<int>{0, 1, 2},
    std::vector<int>{2},
    std::vector<int>{0, 3, 6},
    std::vector<int>{0, 1, 2, 3, 4, 5, 6, 7, 8},
    std::vector<int>{2, 5, 8},
    std::vector<int>{6},
    std::vector<int>{6, 7, 8},
    std::vector<int>{8},

    // second level
    std::vector<int>{0, 9, 18},
    std::vector<int>{0, 1, 2, 9, 10, 11, 18, 19, 20},
    std::vector<int>{2, 11, 20},
    std::vector<int>{0, 3, 6, 9, 12, 15, 18, 21, 24},
    std::vector<int>{0,  1,  2,  3,  4,  5,  6,  7,  8,  9,  10, 11, 12, 13,
                     14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26},
    std::vector<int>{2, 5, 8, 11, 14, 17, 20, 23, 26},
    std::vector<int>{6, 15, 24},
    std::vector<int>{6, 7, 8, 15, 16, 17, 24, 25, 26},
    std::vector<int>{8, 17, 26},

    // third level
    std::vector<int>{18},
    std::vector<int>{18, 19, 20},
    std::vector<int>{20},
    std::vector<int>{18, 21, 24},
    std::vector<int>{18, 19, 20, 21, 22, 23, 24, 25, 26},
    std::vector<int>{20, 23, 26},
    std::vector<int>{24},
    std::vector<int>{24, 25, 26},
    std::vector<int>{26},
};

double eval_integral(double F(double d), double a, double b,
                     double precision = EPSILON) {
  return F(b - precision) - F(a + precision);
}

divVField::divVField() {
  // Integral from -1.5 to -0.5, -0.5 to 0.5, and 0.5 to 1.5 respectively
  dw = {eval_integral(basisF1::div_weight, MIN_OUTER_T, MIN_INNER_T),
        eval_integral(basisF1::div_weight, MIN_INNER_T, MAX_INNER_T),
        eval_integral(basisF1::div_weight, MAX_INNER_T, MAX_OUTER_T)};

  // Integral from -1.5 to -0.5, -0.5 to 0.5, and 0.5 to 1.5 respectively
  wc = {eval_integral(basisF1::div_weight_cmpl, MIN_OUTER_T, MIN_INNER_T),
        eval_integral(basisF1::div_weight_cmpl, MIN_INNER_T, MAX_INNER_T),
        eval_integral(basisF1::div_weight_cmpl, MAX_INNER_T, MAX_OUTER_T)};

  initialize_field(*this);
};

laplaceField::laplaceField() : Field() {
  // Integral from -1.5 to -0.5, -0.5 to 0.5, and 0.5 to 1.5 respectively
  dw = {eval_integral(basisF1::laplace_weight, MIN_OUTER_T, MIN_INNER_T),
        eval_integral(basisF1::laplace_weight, MIN_INNER_T, MAX_INNER_T),
        eval_integral(basisF1::laplace_weight, MAX_INNER_T, MAX_OUTER_T)};

  // Integral from -1.5 to -0.5, -0.5 to 0.5, and 0.5 to 1.5 respectively
  wc = {eval_integral(basisF1::div_weight_cmpl, MIN_OUTER_T, MIN_INNER_T),
        eval_integral(basisF1::div_weight_cmpl, MIN_INNER_T, MAX_INNER_T),
        eval_integral(basisF1::div_weight_cmpl, MAX_INNER_T, MAX_OUTER_T)};

  initialize_field(*this);
};
