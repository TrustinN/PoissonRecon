#include "basis.hpp"
#include <cassert>
#include <functional>

constexpr static double MAX_OUTER_T = 1.5;
constexpr static double MAX_INNER_T = 0.5;
constexpr static double MIN_INNER_T = -0.5;
constexpr static double MIN_OUTER_T = -1.5;
constexpr static double EPSILON = 1e-10;

int sign(double x) { return (x > 0) - (x < 0); }

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

double dbasisF1::operator()(double t) {
  double abs_t = std::abs(t);

  if (abs_t > MAX_OUTER_T) {
    return 0;

  } else if (abs_t > MAX_INNER_T) {
    return -sign(t) * MAX_OUTER_T + abs_t;

  } else {
    return -2 * t;
  };
}

double d2basisF1::operator()(double t) {
  double abs_t = std::abs(t);

  if (abs_t > MAX_OUTER_T) {
    return 0;

  } else if (abs_t > MAX_INNER_T) {
    return 1;

  } else {
    return -2;
  };
}
