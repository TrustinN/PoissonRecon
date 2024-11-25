#ifndef INTEGRATION_HPP
#define INTEGRATION_HPP

#include "utils/linalg.hpp"

#include <array>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_rng.h>

struct node_params {
  std::array<double, 3> c1;
  std::array<double, 3> c2;
  double s;
};

template <typename F0, typename F1, typename F2, typename F3, typename F4,
          typename F5>
double basisFo(double x[], size_t dim, void *params) {
  struct node_params *fp = (struct node_params *)params;

  double t0 = (x[0] - fp->c1[0]) / fp->s;
  double t1 = (x[0] - fp->c2[0]) / fp->s;

  double t2 = (x[1] - fp->c1[1]) / fp->s;
  double t3 = (x[1] - fp->c2[1]) / fp->s;

  double t4 = (x[2] - fp->c1[2]) / fp->s;
  double t5 = (x[2] - fp->c2[2]) / fp->s;

  return F0::operator()(t0) * F1::operator()(t1) * F2::operator()(t2) *
         F3::operator()(t3) * F4::operator()(t4) * F5::operator()(t5);
  ;
}

template <typename F0, typename F1, typename F2, typename F3, typename F4,
          typename F5>
double inner_product(node_params params, double a[], double b[]) {

  gsl_monte_function F;
  gsl_monte_miser_state *workspace = gsl_monte_miser_alloc(3);
  F.f = &basisFo<F0, F1, F2, F3, F4, F5>;
  F.dim = 3;
  F.params = &params;

  double result, error;
  auto T = gsl_rng_default;
  auto r = gsl_rng_alloc(T);

  gsl_monte_miser_integrate(&F, a, b, F.dim, 500000, r, workspace, &result,
                            &error);
  gsl_monte_miser_free(workspace);
  gsl_rng_free(r);

  return result;
}

template <typename F0, typename F1, typename F2, typename F3, typename F4,
          typename F5>
std::array<double, 27> inner_product_field() {

  std::array<double, 3> center{1, 1, 1};
  double s = 0.5;
  std::array<double, 27> res{0};
  std::array<double, 3> l_bounds = {center[0] - 1.5 * s, center[0] - 1.5 * s,
                                    center[0] + .5 * s};
  std::array<double, 3> u_bounds = {center[0] - .5 * s, center[0] + 1.5 * s,
                                    center[0] + 1.5 * s};

  node_params params;
  params.c1 = center;
  params.s = s;
  double length = 2 * s;

  double lower_bound[3];
  double upper_bound[3];

  for (int i = 0; i < 27; i++) {

    std::array<int, 3> axis = {i % 3, (i / 3) % 3, (i / 9) % 3};

    for (int i = 0; i < 3; i++) {
      lower_bound[i] = l_bounds[axis[i]];
      upper_bound[i] = u_bounds[axis[i]];
    }

    std::array<double, 3> f2_center =
        center + std::array<double, 3>{(axis[0] - 1) * length,
                                       (axis[1] - 1) * length,
                                       (axis[2] - 1) * length};

    params.c2 = f2_center;
    res[i] =
        inner_product<F0, F1, F2, F3, F4, F5>(params, lower_bound, upper_bound);
  }
  return res;
};

#endif
