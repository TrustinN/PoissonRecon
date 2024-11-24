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

Field::Field() {
  for (int k = 0; k < 3; k++) {
    for (int j = 0; j < 3; j++) {
      for (int i = 0; i < 3; i++) {
        int idx = i + 3 * j + 9 * k;
        int_field_x[idx] = dw[i] * wc[j] * wc[k];
        int_field_y[idx] = dw[j] * wc[k] * wc[i];
        int_field_z[idx] = dw[k] * wc[i] * wc[j];
      }
    }
  }
};

divVField::divVField() : Field() {};
laplaceField::laplaceField() : Field() {};
