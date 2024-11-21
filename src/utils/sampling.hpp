#ifndef SAMPLING_HPP
#define SAMPLING_HPP

#include <array>
#include <vector>

// -------------------------------------------------------------------------------------------------//
// RANDOM NUMBER GENERATION
// -------------------------------------------------------------------------------------------------//

std::vector<int> rand_ints(int min, int max, int num);
std::vector<std::array<double, 3>> rand_points(double min, double max,
                                               int num_points);

// -------------------------------------------------------------------------------------------------//
// Spatial sampling
// -------------------------------------------------------------------------------------------------//

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
