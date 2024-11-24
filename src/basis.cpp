#include "basis.hpp"
#include <cassert>

double basisF1::div_weight(double t) {
  double abs_t = abs(t);
  if (abs_t > 1.5) {
    return 0;
  } else if (0.5 <= abs_t && abs_t <= 1.5) {
    return std::pow(1.5 - abs_t, 4) / 8;
  } else if (abs_t < 0.5) {
    return -.75 * std::pow(t, 2) + std::pow(t, 4) / 2;
  };
  return 0;
};

double basisF1::div_weight_cmpl(double t) {
  double abs_t = abs(t);
  if (abs_t > 1.5) {
    return 0;
  } else if (0.5 < abs_t && abs_t <= 1.5) {
    return std::pow(1.5 - abs_t, 5) / 20;
  } else if (abs_t <= 0.5) {
    double t_3 = std::pow(t, 3);
    double t_5 = std::pow(t, 5);
    double res = 0.5625 * t - .5 * t_3 + t_5 / 5;
    return res;
  };
  return 0;
};

double basisF1::laplace_weight(double t) {
  double abs_t = abs(t);
  if (abs_t > 1.5) {
    return 0;
  } else if (0.5 < abs_t && abs_t <= 1.5) {
    return 2 * t;
  } else if (abs_t <= 0.5) {
    return -2 * t;
  };
  return 0;
}

double basisF1::operator()(double p) {
  double abs_p = std::abs(p);
  if (abs_p > 1.5) {
    return 0;
  } else if (abs_p > .5) {
    return std::pow(1.5 - abs_p, 2) / 2;
  } else {
    return .75 - std::pow(p, 2);
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
    std::vector<int>{13},
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

divVField::divVField() {

  double d1 = basisF1::div_weight(0.5);
  double d2 = basisF1::div_weight(1.5);

  // Integral from -1.5 to -0.5, -0.5 to 0.5, and 0.5 to 1.5 respectively
  dw = {d1 - d2, 0.0, d2 - d1};

  double dc0 = basisF1::div_weight_cmpl(-0.5);
  double dc1 = basisF1::div_weight_cmpl(0.5);
  double dc2 = basisF1::div_weight_cmpl(1.5);

  // Integral from -1.5 to -0.5, -0.5 to 0.5, and 0.5 to 1.5 respectively
  wc = {dc1 - dc2, dc1 - dc0, dc2 - dc1};
  initialize_field(*this);
};

laplaceField::laplaceField() : Field() {
  double lw1 = basisF1::laplace_weight(-1.5);
  double lw2 = basisF1::laplace_weight(-0.5);
  double lw3 = basisF1::laplace_weight(0.5);
  double lw4 = basisF1::laplace_weight(1.5);

  // Integral from -1.5 to -0.5, -0.5 to 0.5, and 0.5 to 1.5 respectively
  dw = {lw2 - lw1, lw3 - lw2, lw4 - lw3};

  double dc0 = basisF1::div_weight_cmpl(-0.5);
  double dc1 = basisF1::div_weight_cmpl(0.5);
  double dc2 = basisF1::div_weight_cmpl(1.5);

  // Integral from -1.5 to -0.5, -0.5 to 0.5, and 0.5 to 1.5 respectively
  wc = {dc1 - dc2, dc1 - dc0, dc2 - dc1};
  initialize_field(*this);
};
