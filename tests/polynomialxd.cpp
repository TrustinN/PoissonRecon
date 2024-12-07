#include "../src/PolynomialXd.hpp"
#include <gtest/gtest.h>

TEST(PolynomialXdEval, Basic) {
  Polynomial<3> p1({1.0, 2.0, -3, -1});
  Polynomial<3> p2({-2.0, .5, 0, 0});
  Polynomial<3> p3({1.0, 1.5, -3.2, -.2});
  PolynomialXD<3> p({p1, p2, p3});

  std::array<double, 3> q1 = {1, .2, .6};
  double tolerance = 1e-9;
  ASSERT_NEAR(p1(q1[0]) * p2(q1[1]) * p3(q1[2]), p(q1), tolerance);
}
