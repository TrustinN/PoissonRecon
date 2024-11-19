#ifndef NORMAL_HPP
#define NORMAL_HPP

#include <array>
#include <set>
#include <vector>

struct TangentPlane {
  std::array<double, 3> normal;
  std::array<double, 3> center;
};

double offset(const TangentPlane &lhs, const TangentPlane &rhs);

TangentPlane get_tP(std::vector<std::array<double, 3>> vertices);

class NormalApproximations {
private:
  std::vector<std::array<double, 3>> _vertices;
  std::vector<TangentPlane> _planes;
  std::vector<std::set<int>> traversal_order;
  std::vector<std::set<int>> adj_list;

public:
  NormalApproximations(std::vector<std::array<double, 3>> vertices);

  std::vector<std::array<double, 3>> vertices() const { return _vertices; };
  std::vector<TangentPlane> planes() const { return _planes; };
};

#endif
