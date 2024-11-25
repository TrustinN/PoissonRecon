#ifndef LINALG_HPP
#define LINALG_HPP

#include "../Octree.hpp"

// -------------------------------------------------------------------------------------------------//
// LINEAR ALGEBRA
// -------------------------------------------------------------------------------------------------//

std::array<double, 3> operator*(const std::array<double, 3> &a,
                                const double &b);
std::array<double, 3> operator*(const double &a,
                                const std::array<double, 3> &b);
std::array<double, 3> operator+(const std::array<double, 3> &a,
                                const double &b);
std::array<double, 3> operator+(const double &a,
                                const std::array<double, 3> &b);
std::array<double, 3> operator+(const std::array<double, 3> &a,
                                const std::array<double, 3> &b);
std::array<double, 3> operator-(const std::array<double, 3> &a,
                                const std::array<double, 3> &b);
double dot(const std::array<double, 3> &a, const std::array<double, 3> &b);
double dot(const std::array<int, 3> &a, const std::array<double, 3> &b);
double distance(const std::array<double, 3> &a, const Node *node);
double distance(const std::array<double, 3> &a, const std::array<double, 3> &b);

#endif
