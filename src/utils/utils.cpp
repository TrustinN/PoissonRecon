#include "utils.hpp"
#include <algorithm>
#include <array>
#include <random>

// -------------------------------------------------------------------------------------------------//
// LINEAR ALGEBRA
// -------------------------------------------------------------------------------------------------//

std::array<double, 3> operator*(const std::array<double, 3> &a,
                                const double &b) {
  return std::array<double, 3>{a[0] * b, a[1] * b, a[2] * b};
};

std::array<double, 3> operator*(const double &a,
                                const std::array<double, 3> &b) {
  return b * a;
};

std::array<double, 3> operator+(const std::array<double, 3> &a,
                                const double &b) {
  return std::array<double, 3>{a[0] + b, a[1] + b, a[2] + b};
};

std::array<double, 3> operator+(const double &a,
                                const std::array<double, 3> &b) {
  return b + a;
};

std::array<double, 3> operator+(const std::array<double, 3> &a,
                                const std::array<double, 3> &b) {
  return std::array<double, 3>{a[0] + b[0], a[1] + b[1], a[2] + b[2]};
};

std::array<double, 3> operator-(const std::array<double, 3> &a,
                                const std::array<double, 3> &b) {
  return std::array<double, 3>{a[0] - b[0], a[1] - b[1], a[2] - b[2]};
};

double dot(const std::array<double, 3> &a, const std::array<double, 3> &b) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

// helper function computes distance between
double distance(const std::array<double, 3> &a, const Node *node) {
  std::array<double, 3> center = node->center;
  std::array<double, 3> diff = {std::abs(a[0] - center[0]),
                                std::abs(a[1] - center[1]),
                                std::abs(a[2] - center[2])};

  // If distance along an axis is < node width
  // only need to compute remaining distance along other axis
  double width = node->width;
  diff[0] = (diff[0] < width) ? 0 : std::pow(diff[0] - width, 2);
  diff[1] = (diff[1] < width) ? 0 : std::pow(diff[1] - width, 2);
  diff[2] = (diff[2] < width) ? 0 : std::pow(diff[2] - width, 2);

  return diff[0] + diff[1] + diff[2];
};

double distance(const std::array<double, 3> &a,
                const std::array<double, 3> &b) {
  return std::pow(a[0] - b[0], 2) + std::pow(a[1] - b[1], 2) +
         std::pow(a[2] - b[2], 2);
};

// -------------------------------------------------------------------------------------------------//
// GRAPH UTILS
// -------------------------------------------------------------------------------------------------//

std::vector<std::set<int>> join_graphs(const std::vector<std::set<int>> &g1,
                                       const std::vector<std::set<int>> &g2) {
  std::vector<std::set<int>> ret_adj_list(g1.size());
  for (int i = 0; i < g1.size(); i++) {
    std::set_union(g1[i].begin(), g1[i].end(), g2[i].begin(), g2[i].end(),
                   std::inserter(ret_adj_list[i], ret_adj_list[i].begin()));
  }
  return ret_adj_list;
};

// -------------------------------------------------------------------------------------------------//
// RANDOM NUMBER GENERATION
// -------------------------------------------------------------------------------------------------//

std::vector<int> rand_ints(int min, int max, int num) {

  std::random_device rd;
  std::mt19937 gen(rd());

  std::uniform_int_distribution<> dis(min, max);
  std::vector<int> points(num);
  std::generate(points.begin(), points.end(), [&]() { return dis(gen); });
  return points;
};

std::vector<std::array<double, 3>> rand_points(double min, double max,
                                               int num_points) {
  std::random_device rd;
  std::mt19937 gen(rd());

  std::uniform_real_distribution<> dis(min, max);

  std::vector<std::array<double, 3>> points(num_points);
  for (int i = 0; i < num_points; i++) {
    std::generate(points[i].begin(), points[i].end(),
                  [&]() { return dis(gen); });
  };
  return points;
};
