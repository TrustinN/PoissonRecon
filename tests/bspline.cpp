#include "../src/BSpline.hpp"
#include <gtest/gtest.h>

TEST(BSplineIntegral, Basic) {
  double width = 1.0;
  double center = 1.0;
  PPolynomial<2> bs = BSpline;

  double n_center = center - 2 * width;

  ScalarField<2> sf1(bs);
  ScalarField<2> sf2(bs, std::array<double, 3>{n_center, 0, 0});

  // ASSERT_EQ(sf1.innerProduct(sf2), 0.00833333 * std::pow(.55, 2));
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
  std::cout << sf1 << std::endl;
  std::cout << sf2 << std::endl;

  ASSERT_EQ(sf1.innerProduct(sf2), 0.0);
}
