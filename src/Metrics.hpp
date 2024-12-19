#ifndef METRICS_HPP
#define METRICS_HPP

#include "utils/linalg.hpp"
#include <array>
#include <iostream>

struct Metric {
  double threshold = std::numeric_limits<double>::infinity();
  Metric() {};
  Metric(double r) : threshold(r) {};
  virtual bool operator()(const std::array<double, 3> &t1,
                          const std::array<double, 3> &t2) const {
    return true;
  };
};

struct DefaultMetric : public Metric {
  DefaultMetric() {};
  DefaultMetric(double r) : Metric(r) {};

  inline bool operator()(const std::array<double, 3> &t1,
                         const std::array<double, 3> &t2) const override {
    return true;
  }
};

struct WithinRadius : public Metric {
  WithinRadius(double r) : Metric(r) {};

  inline bool operator()(const std::array<double, 3> &t1,
                         const std::array<double, 3> &t2) const override {
    auto diff = t1 - t2;
    return (std::pow(diff[0], 2) + std::pow(diff[1], 2) +
            std::pow(diff[2], 2)) <= std::pow(threshold, 2);
  };
};

struct WithinBoundingBox : public Metric {
  WithinBoundingBox(double r) : Metric(r) {};

  inline bool operator()(const std::array<double, 3> &t1,
                         const std::array<double, 3> &t2) const override {
    auto diff = t1 - t2;
    return (std::abs(diff[0]) <= threshold) &&
           (std::abs(diff[1]) <= threshold) && (std::abs(diff[2]) <= threshold);
  };
};

#endif
