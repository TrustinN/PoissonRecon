#ifndef LINALG_HPP
#define LINALG_HPP

#include <array>

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
double distance(const std::array<double, 3> &a, const std::array<double, 3> &b);

#endif
