#ifndef IO_HPP
#define IO_HPP

#include "../Octree.hpp"
#include "Eigen/Dense"
#include "Eigen/Sparse"
#include <array>
#include <fstream>
#include <iostream>

template <typename T, size_t N>
std::ostream &operator<<(std::ostream &ofs, const std::array<T, N> &a) {
  ofs << "[";
  if (N > 0) {
    ofs << a[0];
    for (int i = 1; i < N; i++) {
      ofs << ", " << a[i];
    }
  };
  return ofs << "]";
};

std::ostream &operator<<(std::ostream &ofs, const id_point &a);

template <typename T>
std::ostream &operator<<(std::ostream &ofs, const std::pair<T, T> &v) {
  return ofs << "(" << std::get<0>(v) << ", " << std::get<1>(v) << ")";
}

template <typename T>
std::ostream &operator<<(std::ostream &ofs, const std::vector<T> &v) {
  ofs << "<";
  if (!v.empty()) {
    ofs << v[0];
    for (std::size_t i = 1; i < v.size(); i++) {
      ofs << ", " << v[i];
    };
  }
  ofs << ">";
  return ofs;
}

std::ostream &operator<<(std::ostream &ofs, const Node &n);
std::ostream &operator<<(std::ostream &ofs, Node *n);
std::ostream &operator<<(std::ostream &ofs, const Octree &o);

void writeVectorToFile(const Eigen::VectorXd &vec, const std::string &filename);
void writeVectorToFile(const std::vector<double> &vec,
                       const std::string &filename);
std::vector<double> loadVectorFromFile(const std::string &filename);

void save_sparse_matrix(const Eigen::SparseMatrix<double> &matrix,
                        const std::string &filename);

void save_points(const std::vector<std::array<double, 3>> &points,
                 const std::string &filename);
std::vector<std::array<double, 3>> load_points(const std::string &filename);

#endif
