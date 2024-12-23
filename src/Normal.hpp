#ifndef NORMAL_HPP
#define NORMAL_HPP

#include <array>
#include <set>
#include <vector>

double offset(const std::array<double, 3> &n1, const std::array<double, 3> &n2);
std::array<double, 3>
get_normal(const std::vector<std::array<double, 3>> &vertices);

// Align n2 to n1
void align_normals(const std::array<double, 3> &n1, std::array<double, 3> &n2);

void dfs_align(int cur_node, const std::vector<std::set<int>> &adj_list,
               const std::vector<int> &visited,
               std::vector<std::array<double, 3>> &normals);

// Modify normals vector to have correctly oriented normals
void orient_normals(std::vector<std::array<double, 3>> &normals,
                    const std::vector<std::array<double, 3>> &centers,
                    std::vector<std::set<int>> traversal_order);

class NormalApproximations {
private:
  std::vector<std::array<double, 3>> _vertices;
  std::vector<std::array<double, 3>> _normals;
  std::vector<std::array<double, 3>> _inward_normals;
  std::vector<std::set<int>> _traversal_order;
  std::vector<std::set<int>> _adj_list;

public:
  NormalApproximations() {};
  NormalApproximations(std::vector<std::array<double, 3>> vertices);

  std::vector<std::array<double, 3>> vertices() const { return _vertices; };
  std::vector<std::array<double, 3>> normals() const { return _normals; };
  std::vector<std::array<double, 3>> inward_normals() const {
    return _inward_normals;
  };
};

#endif
