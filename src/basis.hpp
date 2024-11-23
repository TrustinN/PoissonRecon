#ifndef BASIS_HPP
#define BASIS_HPP

#include <array>

// 1-D basis function
struct basisF1 {

  // gives the integral over t of dF1(t)/dt * F1(t)
  static double div_weight(double t);
  // gives the integral over t of F1^2(t)
  // computed as the complement to the current integral
  static double div_weight_cmpl(double t);
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

// We output the integral of dF1(t_x)/dt_x * F1(t_x) * (F1(t_y) * F1(t_z))^2 wrt
// t_x, t_y, t_z
//
// array is divided into sections of 3 x 3 x 3 = 27
// on intervals [-1.5, -0.5], [-0.5, 0.5], [0.5, 1.5]
struct intVField {
  double d1 = basisF1::div_weight(0.5);
  double d2 = basisF1::div_weight(1.5);

  // Integral from -1.5 to -0.5, -0.5 to 0.5, and 0.5 to 1.5 respectively
  std::array<double, 3> dw = {d1 - d2, 0.0, d2 - d1};

  double dc0 = basisF1::div_weight_cmpl(-0.5);
  double dc1 = basisF1::div_weight_cmpl(0.5);
  double dc2 = basisF1::div_weight_cmpl(1.5);

  // Integral from -1.5 to -0.5, -0.5 to 0.5, and 0.5 to 1.5 respectively
  std::array<double, 3> wc = {dc1 - dc2, dc1 - dc0, dc2 - dc1};

  intVField();
  std::array<double, 27> int_field_x;
  std::array<double, 27> int_field_y;
  std::array<double, 27> int_field_z;
  double _sum_value;
};

#endif
