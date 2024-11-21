#include "sampling.hpp"
#include "utils.hpp"
#include <Eigen/Dense>
#include <iostream>
#include <random>

std::vector<std::array<double, 3>> sample_sphere(int n, double r) {

  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<double> dis_scale(0.0, 1.0);

  std::vector<std::array<double, 3>> points(n);
  for (int i = 0; i < n; i++) {
    double z = 2 * r * dis_scale(gen) - r;
    double theta = 2 * M_PI * dis_scale(gen);

    double a = std::sqrt(std::pow(r, 2) - std::pow(z, 2));
    double x = std::cos(theta) * a;
    double y = std::sin(theta) * a;

    points[i] = {x, y, z};
  };
  return points;
};

std::vector<std::array<double, 3>> sample_tri(int n,
                                              const std::array<double, 3> &c1,
                                              const std::array<double, 3> &c2,
                                              const std::array<double, 3> &c3) {

  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<double> dis_scale(0.0, 1.0);

  std::vector<std::array<double, 3>> points(n);

  for (int i = 0; i < n; i++) {
    double r1 = dis_scale(gen);
    double r2 = dis_scale(gen);

    double lambda1 = 1.0 - sqrt(r1);
    double lambda2 = sqrt(r1) * (1.0 - r2);
    double lambda3 = sqrt(r1) * r2;

    std::array<double, 3> point = {
        lambda1 * c1[0] + lambda2 * c2[0] + lambda3 * c3[0],
        lambda1 * c1[1] + lambda2 * c2[1] + lambda3 * c3[1],
        lambda1 * c1[2] + lambda2 * c2[2] + lambda3 * c3[2]};

    points[i] = point;
  };
  return points;
};

std::vector<std::array<double, 3>>
sample_plane(int n, const std::array<double, 3> &c1,
             const std::array<double, 3> &c2, const std::array<double, 3> &c3,
             const std::array<double, 3> &c4) {

  std::vector<std::array<double, 3>> points = sample_tri(n / 2, c1, c2, c3);
  std::vector<std::array<double, 3>> points2 = sample_tri(n / 2, c1, c3, c4);
  points.insert(points.end(), points2.begin(), points2.end());
  return points;
};

std::vector<std::array<double, 3>> sample_box(int n, double x, double y,
                                              double z) {

  std::array<double, 3> corner1 = {-x, -y, -z};
  std::array<double, 3> corner2 = {x, y, z};
  std::vector<std::array<double, 3>> points;

  for (int i = 0; i < 3; i++) {
    // set {-x, -y, -z} as base corner

    // fix the sign of one variable and cycle the sign of the others
    std::vector<std::array<double, 3>> cc_corners;
    for (int idx1 = 0; idx1 < 3; idx1++) {
      std::array<double, 3> nc = corner1;
      int idx2 = (idx1 + 1) % 3;
      nc[idx1] *= -1;
      nc[idx2] *= -1;
      nc[i] = -std::abs(nc[i]); // fix the sign of the current index to negative
      cc_corners.push_back(nc);
    }
    cc_corners.insert(cc_corners.begin() + i, corner1);
    std::vector<std::array<double, 3>> np = sample_plane(
        n / 8, cc_corners[0], cc_corners[1], cc_corners[2], cc_corners[3]);

    points.insert(points.end(), np.begin(), np.end());

    std::vector<std::array<double, 3>> cc2_corners;
    for (int idx1 = 0; idx1 < 3; idx1++) {
      std::array<double, 3> nc = corner2;
      int idx2 = (idx1 + 1) % 3;
      nc[idx1] *= -1;
      nc[idx2] *= -1;
      nc[i] = std::abs(nc[i]); // fix the sign of the current index to negative
      cc2_corners.push_back(nc);
    }
    cc2_corners.insert(cc2_corners.begin() + i, corner2);
    np = sample_plane(n / 8, cc2_corners[0], cc2_corners[1], cc2_corners[2],
                      cc2_corners[3]);

    points.insert(points.end(), np.begin(), np.end());
  }

  return points;
};
