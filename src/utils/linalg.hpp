#ifndef LINALG_HPP
#define LINALG_HPP

#include "../Octree.hpp"
#include <array>
#include <cmath>

// -------------------------------------------------------------------------------------------------//
// LINEAR ALGEBRA
// -------------------------------------------------------------------------------------------------//

// Addition
template <typename T, std::size_t N>
std::array<T, N> operator+(const std::array<T, N> &a,
                           const std::array<T, N> &b) {
  std::array<T, N> result;
  for (std::size_t i = 0; i < N; ++i) {
    result[i] = a[i] + b[i];
  }
  return result;
}

template <typename T, std::size_t N>
std::array<T, N> &operator+=(std::array<T, N> &a, const std::array<T, N> &b) {
  for (std::size_t i = 0; i < N; ++i) {
    a[i] += b[i];
  }
  return a;
}

// Subtraction
template <typename T, std::size_t N>
std::array<T, N> operator-(const std::array<T, N> &a,
                           const std::array<T, N> &b) {
  std::array<T, N> result;
  for (std::size_t i = 0; i < N; ++i) {
    result[i] = a[i] - b[i];
  }
  return result;
}

template <typename T, std::size_t N>
std::array<T, N> &operator-=(std::array<T, N> &a, const std::array<T, N> &b) {
  for (std::size_t i = 0; i < N; ++i) {
    a[i] -= b[i];
  }
  return a;
}

// Scalar Multiplication
template <typename T, typename Scalar, std::size_t N>
std::array<T, N> operator*(const std::array<T, N> &arr, const Scalar &scalar) {
  std::array<T, N> result;
  for (std::size_t i = 0; i < N; ++i) {
    result[i] = arr[i] * scalar;
  }
  return result;
}

template <typename T, typename Scalar, std::size_t N>
std::array<T, N> &operator*=(std::array<T, N> &arr, const Scalar &scalar) {
  for (std::size_t i = 0; i < N; ++i) {
    arr[i] *= scalar;
  }
  return arr;
}

template <typename T, typename Scalar, std::size_t N>
std::array<T, N> operator*(const Scalar &scalar, const std::array<T, N> &arr) {
  return arr * scalar; // Delegate to the other overload
}

// Scalar Division
template <typename T, typename Scalar, std::size_t N>
std::array<T, N> operator/(const std::array<T, N> &arr, const Scalar &scalar) {
  std::array<T, N> result;
  for (std::size_t i = 0; i < N; ++i) {
    result[i] = arr[i] / scalar;
  }
  return result;
}

template <typename T, typename Scalar, std::size_t N>
std::array<T, N> &operator/=(std::array<T, N> &arr, const Scalar &scalar) {
  for (std::size_t i = 0; i < N; ++i) {
    arr[i] /= scalar;
  }
  return arr;
}

// Dot Product
template <typename T, std::size_t N>
T dot(const std::array<T, N> &a, const std::array<T, N> &b) {
  T result = T{};
  for (std::size_t i = 0; i < N; ++i) {
    result += a[i] * b[i];
  }
  return result;
}

// Magnitude (Euclidean Norm)
template <typename T, std::size_t N> T magnitude(const std::array<T, N> &arr) {
  return std::sqrt(dot(arr, arr));
}

// Normalization
template <typename T, std::size_t N>
std::array<T, N> normalize(const std::array<T, N> &arr) {
  T mag = magnitude(arr);
  if (mag == T{}) {
    throw std::runtime_error("Cannot normalize a zero-length vector.");
  }
  return arr / mag;
}

// Cross Product (only for 3D vectors)
template <typename T>
std::array<T, 3> cross(const std::array<T, 3> &a, const std::array<T, 3> &b) {
  return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2],
          a[0] * b[1] - a[1] * b[0]};
}

double distance(const std::array<double, 3> &a, const Node *node);

template <typename T, std::size_t N>
T distance(const std::array<T, N> &a, const std::array<T, N> &b) {
  T dist = T{};
  for (int i = 0; i < N; ++i) {
    dist += std::pow(a[i] - b[i], 2);
  }
  return dist;
};

#endif
