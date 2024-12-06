#include "../src/PPolynomial.hpp"
#include <gtest/gtest.h>

TEST(PPolynomialConstruct, Default) {
  PPolynomial<3> p;
  Polynomial<3> q;
  for (int i = 0; i < p.intervals.size(); i++) {
    ASSERT_EQ(p.polys[i], q);
  }
}

TEST(PPolynomialConstruct, AddDefault) {
  PPolynomial<3> p;
  Polynomial<3> q1({1});
  Polynomial<3> q2({1, 2});
  Polynomial<3> q3({0, 1, 1, 3});
  std::vector<double> interval{-std::numeric_limits<double>::infinity(), 1, 4};
  PPolynomial<3> q({q1, q2, q3}, interval);
  PPolynomial<3> s = q + p;
  ASSERT_EQ(q, s);
}

TEST(PPolynomialConstruct, NonDefault) {
  Polynomial<3> q1({1});
  Polynomial<3> q2({1, 2});
  Polynomial<3> q3({0, 1, 1, 3});
  std::vector<double> interval{-std::numeric_limits<double>::infinity(), 1, 4};
  PPolynomial<3> p({q1, q2, q3}, interval);
  ASSERT_EQ(p.polys[0], q1);
  ASSERT_EQ(p.polys[1], q2);
  ASSERT_EQ(p.polys[2], q3);
}

TEST(PPolynomialAdd, Zero) {
  Polynomial<3> q1({0, 0, 0, 0});
  Polynomial<3> q2({0, 0, 0, 0});
  Polynomial<3> q3({0, 0, 0, 0});
  std::vector<double> interval{-std::numeric_limits<double>::infinity(), 1, 4};
  PPolynomial<3> q({q1, q2, q3}, interval);

  Polynomial<3> p1({1, 1, 3, 2});
  Polynomial<3> p2({1, 2, 0, 0});
  Polynomial<3> p3({0, 1, 5, 3});
  PPolynomial<3> p({p1, p2, p3}, interval);

  PPolynomial<3> actual_p = p + q;
  std::vector<double> expected_interval = interval;

  Polynomial<3> expected_p1 = p1 + q1;
  Polynomial<3> expected_p2 = p2 + q2;
  Polynomial<3> expected_p3 = p3 + q3;

  std::vector<Polynomial<3>> expected_polys = {expected_p1, expected_p2,
                                               expected_p3};
  PPolynomial<3> expected(expected_polys, expected_interval);

  ASSERT_EQ(actual_p, expected);
}

TEST(PPolynomialAdd, SameInterval) {
  Polynomial<3> q1({1, 0, 1, 0});
  Polynomial<3> q2({1, 2, 3, 4});
  Polynomial<3> q3({0, 1, 1, 3});
  std::vector<double> interval{-std::numeric_limits<double>::infinity(), 1, 4};
  PPolynomial<3> q({q1, q2, q3}, interval);

  Polynomial<3> p1({1, 1, 3, 2});
  Polynomial<3> p2({1, 2, 0, 0});
  Polynomial<3> p3({0, 1, 5, 3});
  PPolynomial<3> p({p1, p2, p3}, interval);

  PPolynomial<3> actual_p = p + q;
  std::vector<double> expected_interval = interval;

  Polynomial<3> expected_p1 = p1 + q1;
  Polynomial<3> expected_p2 = p2 + q2;
  Polynomial<3> expected_p3 = p3 + q3;

  std::vector<Polynomial<3>> expected_polys = {expected_p1, expected_p2,
                                               expected_p3};
  PPolynomial<3> expected(expected_polys, expected_interval);

  ASSERT_EQ(actual_p, expected);
}

TEST(PPolynomialAdd, NonAlignedInterval) {
  Polynomial<3> q1({1, 0, 1, 0});
  Polynomial<3> q2({1, 2, 3, 4});
  Polynomial<3> q3({0, 1, 1, 3});
  std::vector<double> interval1{-std::numeric_limits<double>::infinity(), 1, 4};
  PPolynomial<3> q({q1, q2, q3}, interval1);

  Polynomial<3> p1({1, 1, 3, 2});
  Polynomial<3> p2({1, 2, 0, 0});
  Polynomial<3> p3({0, 1, 5, 3});
  std::vector<double> interval2{-std::numeric_limits<double>::infinity(), 3, 4};
  PPolynomial<3> p({p1, p2, p3}, interval2);

  PPolynomial<3> actual_p = p + q;
  std::vector<double> expected_interval = {
      -std::numeric_limits<double>::infinity(), 1, 3, 4};

  Polynomial<3> expected_p1 = p1 + q1;
  Polynomial<3> expected_p2 = p1 + q2;
  Polynomial<3> expected_p3 = p2 + q2;
  Polynomial<3> expected_p4 = p3 + q3;

  std::vector<Polynomial<3>> expected_polys = {expected_p1, expected_p2,
                                               expected_p3, expected_p4};
  PPolynomial<3> expected(expected_polys, expected_interval);

  ASSERT_EQ(actual_p, expected);
}

TEST(PPolynomialAddConst, Zero) {
  Polynomial<3> p1({1, 1, 3, 2});
  Polynomial<3> p2({1, 2, 0, 0});
  Polynomial<3> p3({0, 1, 5, 3});
  std::vector<double> interval{-std::numeric_limits<double>::infinity(), 1, 4};
  PPolynomial<3> p({p1, p2, p3}, interval);

  double s = 0.0;

  PPolynomial<3> actual_p = p + s;
  ASSERT_EQ(actual_p, p);
  p -= s;
  ASSERT_EQ(actual_p, p);
}

TEST(PPolynomialAddConst, NonZero) {
  Polynomial<3> p1({1, 1, 3, 2});
  Polynomial<3> p2({1, 2, 0, 0});
  Polynomial<3> p3({0, 1, 5, 3});
  std::vector<double> interval{-std::numeric_limits<double>::infinity(), 1, 4};
  PPolynomial<3> p({p1, p2, p3}, interval);

  Polynomial<3> q1({4, 1, 3, 2});
  Polynomial<3> q2({4, 2, 0, 0});
  Polynomial<3> q3({3, 1, 5, 3});
  PPolynomial<3> q({q1, q2, q3}, interval);

  double s = 3.0;

  PPolynomial<3> actual_p = p + s;
  ASSERT_EQ(actual_p, q);
  p += s;
  ASSERT_EQ(p, q);
}

TEST(PPolynomialSub, Zero) {
  Polynomial<3> q1({0, 0, 0, 0});
  Polynomial<3> q2({0, 0, 0, 0});
  Polynomial<3> q3({0, 0, 0, 0});
  std::vector<double> interval{-std::numeric_limits<double>::infinity(), 1, 4};
  PPolynomial<3> q({q1, q2, q3}, interval);

  Polynomial<3> p1({1, 1, 3, 2});
  Polynomial<3> p2({1, 2, 0, 0});
  Polynomial<3> p3({0, 1, 5, 3});
  PPolynomial<3> p({p1, p2, p3}, interval);

  PPolynomial<3> actual_p = p - q;
  std::vector<double> expected_interval = interval;

  Polynomial<3> expected_p1 = p1 - q1;
  Polynomial<3> expected_p2 = p2 - q2;
  Polynomial<3> expected_p3 = p3 - q3;

  std::vector<Polynomial<3>> expected_polys = {expected_p1, expected_p2,
                                               expected_p3};
  PPolynomial<3> expected(expected_polys, expected_interval);

  ASSERT_EQ(actual_p, expected);
}

TEST(PPolynomialSub, SameInterval) {
  Polynomial<3> q1({1, 0, 1, 0});
  Polynomial<3> q2({1, 2, 3, 4});
  Polynomial<3> q3({0, 1, 1, 3});
  std::vector<double> interval{-std::numeric_limits<double>::infinity(), 1, 4};
  PPolynomial<3> q({q1, q2, q3}, interval);

  Polynomial<3> p1({1, 1, 3, 2});
  Polynomial<3> p2({1, 2, 0, 0});
  Polynomial<3> p3({0, 1, 5, 3});
  PPolynomial<3> p({p1, p2, p3}, interval);

  PPolynomial<3> actual_p = p - q;
  std::vector<double> expected_interval = interval;

  Polynomial<3> expected_p1 = p1 - q1;
  Polynomial<3> expected_p2 = p2 - q2;
  Polynomial<3> expected_p3 = p3 - q3;

  std::vector<Polynomial<3>> expected_polys = {expected_p1, expected_p2,
                                               expected_p3};
  PPolynomial<3> expected(expected_polys, expected_interval);

  ASSERT_EQ(actual_p, expected);
}

TEST(PPolynomialSub, NonAlignedInterval) {
  Polynomial<3> q1({1, 0, 1, 0});
  Polynomial<3> q2({1, 2, 3, 4});
  Polynomial<3> q3({0, 1, 1, 3});
  std::vector<double> interval1{-std::numeric_limits<double>::infinity(), 1, 4};
  PPolynomial<3> q({q1, q2, q3}, interval1);

  Polynomial<3> p1({1, 1, 3, 2});
  Polynomial<3> p2({1, 2, 0, 0});
  Polynomial<3> p3({0, 1, 5, 3});
  std::vector<double> interval2{-std::numeric_limits<double>::infinity(), 3, 4};
  PPolynomial<3> p({p1, p2, p3}, interval2);

  PPolynomial<3> actual_p = p - q;
  std::vector<double> expected_interval = {
      -std::numeric_limits<double>::infinity(), 1, 3, 4};

  Polynomial<3> expected_p1 = p1 - q1;
  Polynomial<3> expected_p2 = p1 - q2;
  Polynomial<3> expected_p3 = p2 - q2;
  Polynomial<3> expected_p4 = p3 - q3;

  std::vector<Polynomial<3>> expected_polys = {expected_p1, expected_p2,
                                               expected_p3, expected_p4};
  PPolynomial<3> expected(expected_polys, expected_interval);

  ASSERT_EQ(actual_p, expected);
}

TEST(PPolynomialSubConst, Zero) {
  Polynomial<3> p1({1, 1, 3, 2});
  Polynomial<3> p2({1, 2, 0, 0});
  Polynomial<3> p3({0, 1, 5, 3});
  std::vector<double> interval{-std::numeric_limits<double>::infinity(), 1, 4};
  PPolynomial<3> p({p1, p2, p3}, interval);

  double s = 0.0;

  PPolynomial<3> actual_p = p - s;
  std::vector<double> expected_interval = interval;

  Polynomial<3> expected_p1 = p1 - s;
  Polynomial<3> expected_p2 = p2 - s;
  Polynomial<3> expected_p3 = p3 - s;

  ASSERT_EQ(actual_p, p);
  p -= s;
  ASSERT_EQ(actual_p, p);
}

TEST(PPolynomialSubConst, NonZero) {
  Polynomial<3> p1({1, 1, 3, 2});
  Polynomial<3> p2({1, 2, 0, 0});
  Polynomial<3> p3({0, 1, 5, 3});
  std::vector<double> interval{-std::numeric_limits<double>::infinity(), 1, 4};
  PPolynomial<3> p({p1, p2, p3}, interval);

  Polynomial<3> q1({-2, 1, 3, 2});
  Polynomial<3> q2({-2, 2, 0, 0});
  Polynomial<3> q3({-3, 1, 5, 3});
  PPolynomial<3> q({q1, q2, q3}, interval);

  double s = 3.0;

  PPolynomial<3> actual_p = p - s;
  ASSERT_EQ(actual_p, q);
  p -= s;
  ASSERT_EQ(p, q);
}

TEST(PPolynomialMul, Zero) {
  Polynomial<3> q1({0, 0, 0, 0});
  Polynomial<3> q2({0, 0, 0, 0});
  Polynomial<3> q3({0, 0, 0, 0});
  std::vector<double> interval{-std::numeric_limits<double>::infinity(), 1, 4};
  PPolynomial<3> q({q1, q2, q3}, interval);

  Polynomial<3> p1({1, 1, 3, 2});
  Polynomial<3> p2({1, 2, 0, 0});
  Polynomial<3> p3({0, 1, 5, 3});
  PPolynomial<3> p({p1, p2, p3}, interval);

  PPolynomial<6> actual_p = p * q;
  std::vector<double> expected_interval = interval;

  Polynomial<6> expected_p1 = p1 * q1;
  Polynomial<6> expected_p2 = p2 * q2;
  Polynomial<6> expected_p3 = p3 * q3;

  std::vector<Polynomial<6>> expected_polys = {expected_p1, expected_p2,
                                               expected_p3};
  PPolynomial<6> expected(expected_polys, expected_interval);

  ASSERT_EQ(actual_p, expected);
}

TEST(PPolynomialMul, SameInterval) {
  Polynomial<3> q1({1, 0, 1, 0});
  Polynomial<3> q2({1, 2, 3, 4});
  Polynomial<3> q3({0, 1, 1, 3});
  std::vector<double> interval{-std::numeric_limits<double>::infinity(), 1, 4};
  PPolynomial<3> q({q1, q2, q3}, interval);

  Polynomial<3> p1({1, 1, 3, 2});
  Polynomial<3> p2({1, 2, 0, 0});
  Polynomial<3> p3({0, 1, 5, 3});
  PPolynomial<3> p({p1, p2, p3}, interval);

  PPolynomial<6> actual_p = p * q;
  std::vector<double> expected_interval = interval;

  Polynomial<6> expected_p1 = p1 * q1;
  Polynomial<6> expected_p2 = p2 * q2;
  Polynomial<6> expected_p3 = p3 * q3;

  std::vector<Polynomial<6>> expected_polys = {expected_p1, expected_p2,
                                               expected_p3};
  PPolynomial<6> expected(expected_polys, expected_interval);

  ASSERT_EQ(actual_p, expected);
}

TEST(PPolynomialMul, NonAlignedInterval) {
  Polynomial<3> q1({1, 0, 1, 0});
  Polynomial<3> q2({1, 2, 3, 4});
  Polynomial<3> q3({0, 1, 1, 3});
  std::vector<double> interval1{-std::numeric_limits<double>::infinity(), 1, 4};
  PPolynomial<3> q({q1, q2, q3}, interval1);

  Polynomial<3> p1({1, 1, 3, 2});
  Polynomial<3> p2({1, 2, 0, 0});
  Polynomial<3> p3({0, 1, 5, 3});
  std::vector<double> interval2{-std::numeric_limits<double>::infinity(), 3, 4};
  PPolynomial<3> p({p1, p2, p3}, interval2);

  PPolynomial<6> actual_p = p * q;
  std::vector<double> expected_interval = {
      -std::numeric_limits<double>::infinity(), 1, 3, 4};

  Polynomial<6> expected_p1 = p1 * q1;
  Polynomial<6> expected_p2 = p1 * q2;
  Polynomial<6> expected_p3 = p2 * q2;
  Polynomial<6> expected_p4 = p3 * q3;
  std::vector<Polynomial<6>> expected_polys = {expected_p1, expected_p2,
                                               expected_p3, expected_p4};
  PPolynomial<6> expected(expected_polys, expected_interval);

  ASSERT_EQ(actual_p, expected);
}

TEST(PPolynomialMulConst, Zero) {
  Polynomial<3> p1({1, 1, 3, 2});
  Polynomial<3> p2({1, 2, 0, 0});
  Polynomial<3> p3({0, 1, 5, 3});
  std::vector<double> interval{-std::numeric_limits<double>::infinity(), 1, 4};
  PPolynomial<3> p({p1, p2, p3}, interval);

  Polynomial<3> q1({0, 0, 0, 0});
  Polynomial<3> q2({0, 0, 0, 0});
  Polynomial<3> q3({0, 0, 0, 0});
  PPolynomial<3> q({q1, q2, q3}, interval);

  double t = 0.0;

  PPolynomial<3> actual_p = p * t;
  std::vector<double> expected_interval = interval;

  ASSERT_EQ(actual_p, q);
  p *= t;
  ASSERT_EQ(p, q);
}

TEST(PPolynomialMulConst, NonZero) {
  Polynomial<3> p1({1, 1, 3, 2});
  Polynomial<3> p2({1, 2, 0, 0});
  Polynomial<3> p3({0, 1, 5, 3});
  std::vector<double> interval{-std::numeric_limits<double>::infinity(), 1, 4};
  PPolynomial<3> p({p1, p2, p3}, interval);

  Polynomial<3> q1({3, 3, 9, 6});
  Polynomial<3> q2({3, 6, 0, 0});
  Polynomial<3> q3({0, 3, 15, 9});
  PPolynomial<3> q({q1, q2, q3}, interval);

  double t = 3.0;

  PPolynomial<3> actual_p = p * t;
  std::vector<double> expected_interval = interval;

  ASSERT_EQ(actual_p, q);
  p *= t;
  ASSERT_EQ(p, q);
}

TEST(PPolynomialScale, NonZero) {
  Polynomial<3> p1({1, 1, 3, 2});
  Polynomial<3> p2({1, 2, 0, 0});
  Polynomial<3> p3({0, 1, 5, 3});
  std::vector<double> interval{-std::numeric_limits<double>::infinity(), 1, 4};
  PPolynomial<3> p({p1, p2, p3}, interval);

  double s = 3.0;
  PPolynomial<3> q = p.scale(s);

  for (int i = 0; i < q.polys.size(); i++) {
    ASSERT_EQ(q.polys[i], p.polys[i].scale(s));
  }
}

TEST(PPolynomialShift, NonZero) {
  Polynomial<3> p1({1, 1, 3, 2});
  Polynomial<3> p2({1, 2, 0, 0});
  Polynomial<3> p3({0, 1, 5, 3});
  std::vector<double> interval{-std::numeric_limits<double>::infinity(), 1, 4};
  PPolynomial<3> p({p1, p2, p3}, interval);

  double t = 3.0;
  PPolynomial<3> q = p.shift(t);

  for (int i = 0; i < q.polys.size(); i++) {
    ASSERT_EQ(q.polys[i], p.polys[i].shift(t));
  }
}

TEST(PPolynomialEval, Basic) {
  Polynomial<3> p1({1, 1, 3, 2});
  Polynomial<3> p2({1, 2, 0, 0});
  Polynomial<3> p3({0, 1, 5, 3});
  std::vector<double> interval{-std::numeric_limits<double>::infinity(), 1, 4};
  PPolynomial<3> p({p1, p2, p3}, interval);

  std::array<double, 4> x = {-4, -.7, 1.5, 20};
  std::array<double, 4> expected_y = {p1(x[0]), p1(x[1]), p2(x[2]), p3(x[3])};

  for (int i = 0; i < x.size(); i++) {
    ASSERT_EQ(expected_y[i], p(x[i]));
  }
}

TEST(PPolynomialIntegral, Basic) {
  Polynomial<3> p1({1, 1, 3, 2});
  Polynomial<3> p2({1, 2, 0, 0});
  Polynomial<3> p3({0, 1, 5, 3});
  std::vector<double> interval{-std::numeric_limits<double>::infinity(), 1, 4};
  PPolynomial<3> p({p1, p2, p3}, interval);
  PPolynomial<4> q({p1.integral(), p2.integral(), p3.integral()}, interval);

  ASSERT_EQ(p.integral(), q);
}

TEST(PPolynomialDerivative, Basic) {
  Polynomial<3> p1({1, 1, 3, 2});
  Polynomial<3> p2({1, 2, 0, 0});
  Polynomial<3> p3({0, 1, 5, 3});
  std::vector<double> interval{-std::numeric_limits<double>::infinity(), 1, 4};
  PPolynomial<3> p({p1, p2, p3}, interval);
  PPolynomial<2> q({p1.derivative(), p2.derivative(), p3.derivative()},
                   interval);

  ASSERT_EQ(p.derivative(), q);
}

TEST(PPolynomialIntegralEval, Basic) {
  Polynomial<3> p1({1, 1, 3, 2});
  Polynomial<3> p2({1, 2, 0, 0});
  Polynomial<3> p3({0, 1, 5, 3});
  std::vector<double> interval{-std::numeric_limits<double>::infinity(), 1, 4};
  PPolynomial<3> p({p1, p2, p3}, interval);

  double a = -20;
  double b = -8;
  double actual = p.integral(a, b);
  double expected = p1.integral(a, b);

  ASSERT_EQ(expected, actual);
}

TEST(PPolynomialIntegralEval, MultiInterval) {
  Polynomial<3> p1({1, 1, 3, 2});
  Polynomial<3> p2({1, 2, 0, 0});
  Polynomial<3> p3({0, 1, 5, 3});
  std::vector<double> interval{-std::numeric_limits<double>::infinity(), 1, 4};
  PPolynomial<3> p({p1, p2, p3}, interval);

  double a = -20;
  double b = 3;
  double actual = p.integral(a, b);
  double expected = p1.integral(a, 1) + p2.integral(1, b);

  ASSERT_EQ(expected, actual);
}

TEST(PPolynomialIntegralEval, Edge) {
  Polynomial<3> p1({1, 1, 3, 2});
  Polynomial<3> p2({1, 2, 0, 0});
  Polynomial<3> p3({0, 1, 5, 3});
  Polynomial<3> p4({7, 1, 1, 3});
  std::vector<double> interval{-std::numeric_limits<double>::infinity(), 1, 4,
                               5};
  PPolynomial<3> p({p1, p2, p3, p4}, interval);

  double a = -2;
  double b = 5;
  double actual = p.integral(a, b);
  double expected = p1.integral(a, 1) + p2.integral(1, 4) + p3.integral(4, b);

  ASSERT_EQ(expected, actual);
}

TEST(PPolynomialIntegralEval, FullInterval) {
  Polynomial<3> p1({0, 0, 0, 0});
  Polynomial<3> p2({1, 2, 0, 0});
  Polynomial<3> p3({0, 1, 5, 3});
  Polynomial<3> p4({0, 0, 0, 0});
  std::vector<double> interval{-std::numeric_limits<double>::infinity(), 1, 4,
                               5};
  PPolynomial<3> p({p1, p2, p3, p4}, interval);

  double a = -std::numeric_limits<double>::infinity();
  double b = std::numeric_limits<double>::infinity();
  double actual = p.integral(a, b);
  double expected = p1.integral(a, 1) + p2.integral(1, 4) + p3.integral(4, 5) +
                    p4.integral(5, b);

  ASSERT_EQ(expected, actual);
}
