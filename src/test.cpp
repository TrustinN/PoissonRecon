#include "Metrics.hpp"
#include "utils/linalg.hpp"

bool useMetric(const std::array<double, 3> &c1, const std::array<double, 3> &c2,
               const Metric &metric = DefaultMetric()) {
  return metric(c1, c2);
}

int main() {
  std::array<double, 3> c = {0, 0, 0};
  std::array<double, 3> c2 = {0.1, 0.4, -0.3};
  double r = 0.5;
  WithinBoundingBox metric = WithinBoundingBox(r);
  useMetric(c, c2, metric);

  return 0;
}
