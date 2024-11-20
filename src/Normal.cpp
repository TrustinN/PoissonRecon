#include "Normal.hpp"
#include "Emst.hpp"
#include "Octree.hpp"
#include "RiemannianGraph.hpp"
#include "utils.hpp"
#include <Eigen/Dense>
#include <numeric>

double offset(const std::array<double, 3> &n1,
              const std::array<double, 3> &n2) {
  return 1 - std::abs(dot(n1, n2));
};

TangentPlane get_tp(const std::vector<std::array<double, 3>> &vertices) {
  TangentPlane tp;

  std::array<double, 3> centroid = std::accumulate(
      vertices.begin(), vertices.end(), std::array<double, 3>{0, 0, 0});

  for (auto &val : tp.center) {
    val /= vertices.size();
  }

  Eigen::MatrixXd mat(vertices.size(), 3);
  for (int i = 0; i < vertices.size(); i++) {
    mat.row(i) = Eigen::RowVectorXd::Map(centroid.data(), 3) -
                 Eigen::RowVectorXd::Map(vertices[i].data(), 3);
  };

  Eigen::BDCSVD<Eigen::MatrixXd> svd(mat, Eigen::ComputeThinV);
  Eigen::Matrix<double, Eigen::Dynamic, 3> eig_vec = svd.matrixV().col(2);

  std::copy(eig_vec.data(), eig_vec.data() + 3, tp.normal.begin());
  return tp;
};

void align_normals(const std::array<double, 3> &n1, std::array<double, 3> &n2) {
  if (dot(n1, n2) < 0) {
    n2[0] *= -1;
    n2[1] *= -1;
    n2[2] *= -1;
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

  visited[start_idx] = 1;
  dfs_align(start_idx, traversal_order, visited, normals);
};

NormalApproximations::NormalApproximations(
    std::vector<std::array<double, 3>> vertices) {
  Octree octree(vertices);

  RiemannianGraph rg = RiemannianGraph(vertices, octree, 15);
  std::vector<std::set<int>> rg_adj_list = rg.adj_list();

  // Get the tangent planes
  for (int i = 0; i < vertices.size(); i++) {
    std::set<int> &ns = rg_adj_list[i];
    std::vector<std::array<double, 3>> np;
    for (int n : ns) {
      np.push_back(vertices[n]);
    };

    TangentPlane tp = get_tp(np);
    _normals.push_back(tp.normal);
    _planes.push_back(tp);
  };

  // Get the emst
  Emst emst = Emst(_vertices, octree);

  _adj_list = join_graphs(emst.adj_list(), rg_adj_list);
  _traversal_order =
      get_mst<std::array<double, 3>>(_normals, _adj_list, offset);

  orient_normals(_normals, vertices, _traversal_order);
}
