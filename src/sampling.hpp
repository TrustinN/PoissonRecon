
#ifndef SAMPLING_HPP
#define SAMPLING_HPP

#include <array>
#include <vector>

std::vector<std::array<double, 3>> sample_sphere(int n, double r);
std::vector<std::array<double, 3>> sample_tri(int n,
                                              const std::array<double, 3> &c1,
                                              const std::array<double, 3> &c2,
                                              const std::array<double, 3> &c3);

std::vector<std::array<double, 3>>
sample_plane(int n, const std::array<double, 3> &c1,
             const std::array<double, 3> &c2, const std::array<double, 3> &c3,
             const std::array<double, 3> &c4);

std::vector<std::array<double, 3>> sample_box(int n, double x, double y,
                                              double z);

#endif
