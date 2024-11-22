#include "basis.hpp"

double basisF1::operator()(double p) {
  double abs_p = std::abs(p);
  if (abs_p >= 1.5) {
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
