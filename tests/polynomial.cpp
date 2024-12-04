#include "../src/Polynomial.hpp"
#include "../src/utils/linalg.hpp"
#include <gtest/gtest.h>

TEST(PolynomialConstruct, Empty) {
  Polynomial<3> p;
  auto zeros = std::array<double, 4>();
  ASSERT_EQ(p.coefficients, zeros);
}

TEST(PolynomialConstruct, NonZero) {
  std::array<double, 4> coeff = {1, 2, 3, 4};
  Polynomial<3> p(coeff);
  ASSERT_EQ(p.coefficients, coeff);
}

TEST(PolynomialEqual, Default) {
  Polynomial<3> p1;
  Polynomial<3> p2;
  ASSERT_EQ(p1, p2);
}

TEST(PolynomialEqual, NonZero) {
  std::array<double, 4> coeff = {1, 2, 3, 4};
  Polynomial<3> p1(coeff);
  Polynomial<3> p2(coeff);
  ASSERT_EQ(p1, p2);
}

TEST(PolynomialNotEqual, Size) {
  Polynomial<3> p1;
  Polynomial<4> p2;
  ASSERT_NE(p1, p2);
}

TEST(PolynomialNotEqual, Coefficients) {
  std::array<double, 4> coeff1 = {1, 2, 3, 4};
  Polynomial<3> p1(coeff1);
  std::array<double, 4> coeff2 = {1, 1, 3, 4};
  Polynomial<3> p2(coeff2);
  ASSERT_NE(p1, p2);
}

TEST(PolynomialAdd, Zero) {
  std::array<double, 4> coeff1 = {1, 2, 3, 4};
  Polynomial<3> expected(coeff1);

  std::array<double, 4> coeff2 = {0};
  Polynomial<3> p2(coeff2);

  auto actual = expected + p2;

  ASSERT_EQ(actual, expected);
}

TEST(PolynomialAddPoly, Basic) {
  std::array<double, 4> coeff1 = {1, 2, 3, 4};
  Polynomial<3> p1(coeff1);

  std::array<double, 4> coeff2 = {1, 1, 1, 2};
  Polynomial<3> p2(coeff2);

  auto actual = p1 + p2;
  auto coeff3 = coeff1 + coeff2;
  Polynomial<3> expected(coeff3);

  ASSERT_EQ(actual, expected);

  p1 += p2;
  ASSERT_EQ(p1, expected);
}

TEST(PolynomialAddConst, Basic) {
  std::array<double, 4> coeff1 = {1, 2, 3, 4};
  Polynomial<3> p1(coeff1);

  double v{4};
  Polynomial<3> expected({5, 2, 3, 4});

  auto actual = p1 + v;

  ASSERT_EQ(actual, expected);
}

TEST(PolynomialSubtract, Basic) {
  std::array<double, 4> coeff1 = {1, 2, 3, 4};
  Polynomial<3> p1(coeff1);

  std::array<double, 4> coeff2 = {1, 1, 1, 2};
  Polynomial<3> p2(coeff2);

  auto actual = p1 + p2;
  auto coeff3 = coeff1 + coeff2;
  Polynomial<3> expected(coeff3);

  ASSERT_EQ(actual, expected);
}
