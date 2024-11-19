#include "Normal.hpp"
#include "Octree.hpp"
#include "RiemannianGraph.hpp"
#include "utils.hpp"

double offset(const TangentPlane &lhs, const TangentPlane &rhs) {
  return 1 - std::abs(dot(lhs.normal, rhs.normal));
};

NormalApproximations::NormalApproximations(
    std::vector<std::array<double, 3>> vertices) {
  Octree octree(vertices);

  RiemannianGraph rg = RiemannianGraph(vertices, octree, 15);
};
