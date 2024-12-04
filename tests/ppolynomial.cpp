#include "../src/PPolynomial.hpp"
#include <gtest/gtest.h>

TEST(PPolynomialConstruct, Default) {
  PPolynomial<3> p;
  Polynomial<3> q;
  for (int i = 0; i < p.intervals.size(); i++) {
    ASSERT_EQ(p.polys[i], q);
  }
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

  ASSERT_EQ(actual_p.intervals, expected_interval);
  ASSERT_EQ(actual_p.polys[0], expected_p1);
  ASSERT_EQ(actual_p.polys[1], expected_p2);
  ASSERT_EQ(actual_p.polys[2], expected_p3);
  ASSERT_EQ(actual_p.polys.size(), 3);
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

  ASSERT_EQ(actual_p.intervals, expected_interval);
  ASSERT_EQ(actual_p.polys[0], expected_p1);
  ASSERT_EQ(actual_p.polys[1], expected_p2);
  ASSERT_EQ(actual_p.polys[2], expected_p3);
  ASSERT_EQ(actual_p.polys.size(), 3);
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

  ASSERT_EQ(actual_p.intervals, expected_interval);
  ASSERT_EQ(actual_p.polys[0], expected_p1);
  ASSERT_EQ(actual_p.polys[1], expected_p2);
  ASSERT_EQ(actual_p.polys[2], expected_p3);
  ASSERT_EQ(actual_p.polys[3], expected_p4);
  ASSERT_EQ(actual_p.polys.size(), 4);
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

  ASSERT_EQ(actual_p.intervals, expected_interval);
  ASSERT_EQ(actual_p.polys[0], expected_p1);
  ASSERT_EQ(actual_p.polys[1], expected_p2);
  ASSERT_EQ(actual_p.polys[2], expected_p3);
  ASSERT_EQ(actual_p.polys.size(), 3);
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

  ASSERT_EQ(actual_p.intervals, expected_interval);
  ASSERT_EQ(actual_p.polys[0], expected_p1);
  ASSERT_EQ(actual_p.polys[1], expected_p2);
  ASSERT_EQ(actual_p.polys[2], expected_p3);
  ASSERT_EQ(actual_p.polys.size(), 3);
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

  ASSERT_EQ(actual_p.intervals, expected_interval);
  ASSERT_EQ(actual_p.polys[0], expected_p1);
  ASSERT_EQ(actual_p.polys[1], expected_p2);
  ASSERT_EQ(actual_p.polys[2], expected_p3);
  ASSERT_EQ(actual_p.polys[3], expected_p4);
  ASSERT_EQ(actual_p.polys.size(), 4);
}
