#include "../src/BSpline.hpp"
#include <gtest/gtest.h>

TEST(BSplineIntegral, Basic) {
  double width = 1.0;
  double center = 1.0;
  PPolynomial<2> bs = BSpline;

  double n_center = center - 2 * width;

  ScalarField<2> sf1(bs);
  ScalarField<2> sf2(bs, std::array<double, 3>{n_center, 0, 0});

  ASSERT_EQ(sf1.innerProduct(sf2), 0.00833333 * std::pow(.55, 2));
}

TEST(BSplineConstruct, Shift) {
  // ----|--------|--------|----|
  // 1.5 1.0      0        1.0  1.5

  std::array<double, 3> center = {2.25, 0, 0};
  ScalarField<2> sf2(BSpline, center);

  ASSERT_EQ(sf2(center), std::pow(0.75, 3));
}

TEST(BSplineConstruct, Scale) {
  //    |--|----|----|--|
  // -.75 -.5   0   .5 .75

  double width = .5;
  std::array<double, 3> center = {0, 0, 0};
  ScalarField<2> sf2(BSpline, center, 1 / width);

  ASSERT_EQ(sf2({0.75, 0.75, 0.75}), 0);
}

TEST(BSplineIntegral, NonOverlappingSupport) {
  // ----|--------|--------|----|
  // 1.5 1.0      0        1.0  1.5
  //
  //                            |--|----|----|--|
  //                          1.5 1.75 2.25 2.75 3

  double width1 = 1.0;
  double center1 = 0.0;
  PPolynomial<2> bs = BSpline;

  ScalarField<2> sf1(bs);

  double width2 = .5;
  std::array<double, 3> center2 = {2.25, 0, 0};
  ScalarField<2> sf2(bs, center2, 1 / width2);

  ASSERT_EQ(sf1.innerProduct(sf2), 0.0);
}

TEST(BSplineIntegral, Overlapping) {
  // ----|--------|--------|----|
  // 1.5 1.0      0        1.0  1.5
  //
  //               |--|----|----|--|
  //              .25 0.5 1.0  1.5 1.75

  double width1 = 1.0;
  double center1 = 0.0;
  PPolynomial<2> bs1 = BSpline;

  double width2 = .5;
  double center2 = 1.0;
  PPolynomial<2> bs2 = basisFFactory(BSpline, center2, 1 / width2);

  bs2 = bs2.derivative_keep_dim();
  PPolynomial<4> bs3 = bs1 * bs2;

  ASSERT_NEAR(bs3.integral(-5, 5), 0.248697916667, 1e-5);
}
