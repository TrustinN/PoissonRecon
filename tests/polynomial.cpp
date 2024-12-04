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

TEST(PolynomialConstruct, Underfilled) {
  Polynomial<3> p({1});
  std::array<double, 4> expected{1, 0, 0, 0};
  ASSERT_EQ(p.coefficients, expected);
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

TEST(PolynomialSub, Zero) {
  std::array<double, 4> coeff1 = {1, 2, 3, 4};
  Polynomial<3> expected(coeff1);

  std::array<double, 4> coeff2 = {0};
  Polynomial<3> p2(coeff2);

  auto actual = expected - p2;

  ASSERT_EQ(actual, expected);
}

TEST(PolynomialSubPoly, Basic) {
  std::array<double, 4> coeff1 = {1, 2, 3, 4};
  Polynomial<3> p1(coeff1);

  std::array<double, 4> coeff2 = {1, 1, 1, 2};
  Polynomial<3> p2(coeff2);

  auto actual = p1 - p2;
  auto coeff3 = coeff1 - coeff2;
  Polynomial<3> expected(coeff3);

  ASSERT_EQ(actual, expected);

  p1 -= p2;
  ASSERT_EQ(p1, expected);
}

TEST(PolynomialSubConst, Basic) {
  std::array<double, 4> coeff1 = {1, 2, 3, 4};
  Polynomial<3> p1(coeff1);

  double v{4};
  Polynomial<3> expected({-3, 2, 3, 4});

  auto actual = p1 - v;

  ASSERT_EQ(actual, expected);
}

TEST(PolynomialMul, Zero) {
  std::array<double, 4> coeff1 = {1, 2, 3, 4};
  Polynomial<3> p1(coeff1);

  std::array<double, 1> coeff2 = {0};
  Polynomial<0> p2(coeff2);

  Polynomial<3> expected({0, 0, 0, 0});
  auto actual = p1 * p2;

  ASSERT_EQ(actual, expected);
}

TEST(PolynomialMulPoly, Basic) {
  std::array<double, 4> coeff1 = {1, 1, 1, 1};
  Polynomial<3> p1(coeff1);

  std::array<double, 2> coeff2 = {-1, 1};
  Polynomial<1> p2(coeff2);

  Polynomial<4> actual = p1 * p2;
  Polynomial<4> expected({-1, 0, 0, 0, 1});

  ASSERT_EQ(actual, expected);
}

TEST(PolynomialMulConst, Basic) {
  std::array<double, 4> coeff1 = {1, 2, 3, 4};
  Polynomial<3> p1(coeff1);

  double v{4};
  Polynomial<3> expected({4, 8, 12, 16});

  auto actual = p1 * v;

  ASSERT_EQ(actual, expected);
}

TEST(PolynomialDiv, Zero) {
  std::array<double, 4> coeff1 = {1, 2, 3, 4};
  Polynomial<3> expected(coeff1);

  std::array<double, 4> coeff2 = {0};
  Polynomial<3> p2(coeff2);

  auto actual = expected - p2;

  ASSERT_EQ(actual, expected);
}

TEST(PolynomialDivConst, Basic) {
  std::array<double, 4> coeff1 = {1, 2, 3, 4};
  Polynomial<3> p1(coeff1);

  double v{4};
  Polynomial<3> expected({1.0 / 4, 1.0 / 2, 3.0 / 4, 1.0});

  auto actual = p1 / v;

  ASSERT_EQ(actual, expected);

  p1 /= v;
  ASSERT_EQ(p1, expected);
}

TEST(PolynomialShift, Zero) {
  std::array<double, 4> coeff1 = {1, 2, 3, 4};
  Polynomial<3> expected(coeff1);

  double shift = 0;
  auto actual = expected.shift(shift);

  ASSERT_EQ(actual, expected);
}

TEST(PolynomialShift, NonZero) {
  std::array<double, 4> coeff1 = {1, 2, 3, 4};
  Polynomial<3> p1(coeff1);

  double shift = 3;
  auto actual = p1.shift(shift);

  Polynomial<3> expected({-86, 92, -33, 4});
  ASSERT_EQ(actual, expected);
}

TEST(PolynomialScale, One) {
  std::array<double, 4> coeff1 = {1, 2, 3, 4};
  Polynomial<3> expected(coeff1);

  double scale = 1;
  auto actual = expected.scale(scale);

  ASSERT_EQ(actual, expected);
}

TEST(PolynomialScale, NonOne) {
  std::array<double, 4> coeff1 = {1, 2, 3, 4};
  Polynomial<3> p1(coeff1);

  double scale = 3;
  auto actual = p1.scale(scale);

  Polynomial<3> expected({1.0, 2.0 / 3.0, 1.0 / 3.0, 4.0 / 27.0});

  ASSERT_EQ(actual, expected);
}

TEST(PolynomialEval, Basic) {
  std::array<double, 4> coeff1 = {1, 2, 3, 4};
  Polynomial<3> p1(coeff1);

  for (int i = 0; i < 4; i++) {
    ASSERT_EQ(p1(i), 1 + 2 * i + 3 * std::pow(i, 2) + 4 * std::pow(i, 3));
  }
}

TEST(PolynomialIntegral, Basic) {
  std::array<double, 4> coeff1 = {1, 2, 3, 4};
  Polynomial<3> p1(coeff1);

  Polynomial<4> actual = p1.integral();
  Polynomial<4> expected({0, 1, 1, 1, 1});
  ASSERT_EQ(actual, expected);
}

TEST(PolynomialDerivative, Basic) {
  std::array<double, 4> coeff1 = {1, 2, 3, 4};
  Polynomial<3> p1(coeff1);

  Polynomial<2> actual = p1.derivative();
  Polynomial<2> expected({2, 6, 12});
  ASSERT_EQ(actual, expected);
}

TEST(PolynomialIntegralEval, Basic) {
  std::array<double, 4> coeff1 = {1, 2, 3, 4};
  Polynomial<3> p1(coeff1);

  double a = 3, b = 5;
  double actual = p1.integral(a, b);
  Polynomial<4> int_p1 = p1.integral();

  double b_i = int_p1(b);
  double a_i = int_p1(a);
  double expected = b_i - a_i;
  ASSERT_EQ(actual, expected);
}
