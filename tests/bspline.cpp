#include "../src/BSpline.hpp"
#include <gtest/gtest.h>

TEST(BSplineIntegral, Basic) {
  double width = 1.0;
  double center = 1.0;
  PPolynomial<2> bs = BSpline;

  double n_center = center - 2 * width;

  ScalarField<2> sf1(bs);
  ScalarField<2> sf2(bs, std::array<double, 3>{n_center, 0, 0});

  std::cout << bs << std::endl;

  ASSERT_EQ(sf1.innerProduct(sf2), 0.00833333 * std::pow(.55, 2));
}
