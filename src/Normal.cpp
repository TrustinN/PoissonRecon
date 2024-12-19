#include "Normal.hpp"
#include "Emst.hpp"
#include "Octree.hpp"
#include "RiemannianGraph.hpp"
#include "utils/graphs.hpp"
#include "utils/linalg.hpp"
#include <Eigen/Dense>

double offset(const std::array<double, 3> &n1,
              const std::array<double, 3> &n2) {
  return 1 - std::abs(dot(n1, n2));
};

std::array<double, 3>
get_normal(const std::vector<std::array<double, 3>> &vertices) {
  std::array<double, 3> centroid = {0, 0, 0};
  for (const auto &v : vertices) {
    centroid += v;
  };

  for (auto &val : centroid) {
    val /= vertices.size();
  }

  Eigen::MatrixXd mat(vertices.size(), 3);
  for (int i = 0; i < vertices.size(); i++) {
    mat.row(i) = Eigen::RowVectorXd::Map(centroid.data(), 3) -
                 Eigen::RowVectorXd::Map(vertices[i].data(), 3);
  };

  Eigen::BDCSVD<Eigen::MatrixXd> svd(mat, Eigen::ComputeThinV);
  Eigen::Vector3d eig_vec = svd.matrixV().col(2);

  return {eig_vec(0), eig_vec(1), eig_vec(2)};
};

inline void align_normals(const std::array<double, 3> &n1,
                          std::array<double, 3> &n2) {
  if (dot(n1, n2) < 0) {
    n2 *= -1;
  };
};

void dfs_align(int cur_node, const std::vector<std::set<int>> &adj_list,
               std::vector<int> &visited,
               std::vector<std::array<double, 3>> &normals) {
  std::vector<int> stack = {cur_node};
  while (!stack.empty()) {
    int v = stack.back();
    stack.pop_back();
    for (int n : adj_list[v]) {
      if (!visited[n]) {
        visited[n] = 1;
        align_normals(normals[v], normals[n]);
        stack.push_back(n);
      }
    }
  }
};

void orient_normals(std::vector<std::array<double, 3>> &normals,
                    const std::vector<std::array<double, 3>> &centers,
                    std::vector<std::set<int>> traversal_order) {
  std::vector<int> visited(normals.size(), 0);
  // Get highest point along one dimension
  double min_x = std::numeric_limits<double>::infinity();
  int start_idx = 0;
  for (int i = 0; i < normals.size(); i++) {
    if (centers[i][0] < min_x) {
      min_x = centers[i][0];
      start_idx = i;
    }
  };
  align_normals(std::array<double, 3>{-1, 0, 0}, normals[start_idx]);

  visited[start_idx] = 1;
  dfs_align(start_idx, traversal_order, visited, normals);
};

NormalApproximations::NormalApproximations(
    std::vector<std::array<double, 3>> vertices) {
  int num_vertices = vertices.size();
  _vertices = vertices;
  _normals = std::vector<std::array<double, 3>>(num_vertices);
  _inward_normals = std::vector<std::array<double, 3>>(num_vertices);

  Octree octree(vertices);
  RiemannianGraph rg = RiemannianGraph(vertices, octree, 15);
  std::vector<std::set<int>> rg_adj_list = rg.adj_list();

  // Get the tangent planes
#pragma omp parallel for
  for (int i = 0; i < num_vertices; i++) {
    std::set<int> &ns = rg_adj_list[i];
    std::vector<std::array<double, 3>> np;
    for (int n : ns) {
      np.push_back(vertices[n]);
    };

    _normals[i] = (get_normal(np));
  };

  // Get the emst
  Emst emst = Emst(vertices, octree);

  _adj_list = join_graphs(emst.adj_list(), rg_adj_list);
  _traversal_order =
      get_mst<std::array<double, 3>>(_normals, _adj_list, offset);

  orient_normals(_normals, vertices, _traversal_order);
#pragma omp parallel for
  for (int i = 0; i < _normals.size(); i++) {
    _inward_normals[i] = _normals[i] * -1;
  }
}
