#include "io.hpp"
#include "../Octree.hpp"
#include <array>
#include <iostream>

std::ostream &operator<<(std::ostream &ofs, const id_point &a) {
  return ofs << std::get<1>(a);
};

std::ostream &operator<<(std::ostream &ofs, const std::array<Node *, 8> &a) {
  if (a[0] == nullptr) {
    ofs << "[null";
  } else {
    ofs << std::endl;
    ofs << *a[0];
  }
  for (int i = 1; i < 8; i++) {
    if (a[i] == nullptr) {
      ofs << ", null";
    } else {
      ofs << std::endl << *a[i];
    }
  };
  if (a[0] == nullptr) {
    ofs << "]";
  }
  return ofs;
};

std::ostream &operator<<(std::ostream &ofs, const Node &n) {
  std::string indent(2 * n.depth, ' ');
  std::string label = (n.is_leaf ? "Leaf" : "Branch");
  ofs << indent << label << "[" << n.num_points << "] {" << std::endl;
  ofs << indent << "  " << "center: " << n.center << "," << std::endl;
  ofs << indent << "  " << "width: " << n.width << "," << std::endl;
  if (n.is_leaf) {
    ofs << indent << "  " << n.children.points << std::endl;
  } else {
    ofs << indent << "  " << "children:" << n.children.nodes << std::endl;
  }
  ofs << indent << "}";
  return ofs;
};

std::ostream &operator<<(std::ostream &ofs, Node *n) {
  std::string indent(2 * n->depth, ' ');
  std::string label = (n->is_leaf ? "Leaf" : "Branch");
  ofs << indent << label << "[" << n->num_points << "] {" << std::endl;
  ofs << indent << "  " << "center: " << n->center << "," << std::endl;
  ofs << indent << "  " << "width: " << n->width << "," << std::endl;
  if (n->is_leaf) {
    ofs << indent << "  " << n->children.points << std::endl;
  } else {
    ofs << indent << "  " << "children:" << n->children.nodes << std::endl;
  }
  ofs << indent << "}";
  return ofs;
};

std::ostream &operator<<(std::ostream &ofs, const Octree &o) {
  ofs << "Octree: [" << std::endl;
  ofs << *o.root() << std::endl;
  ofs << "]";
  return ofs;
}

void writeVectorToFile(const Eigen::VectorXd &vec,
                       const std::string &filename) {
  std::ofstream file(filename);
  if (!file.is_open()) {
    std::cerr << "Failed to open file: " << filename << std::endl;
    return;
  }

  if (vec.size() > 0) {
    file << vec[0];
    for (int i = 1; i < vec.size(); ++i) {
      file << ", " << vec[i];
    }
  }

  file.close();
  if (file.good()) {
    std::cout << "Vector successfully written to " << filename << std::endl;
  } else {
    std::cerr << "Error occurred while writing to the file." << std::endl;
  }
}

void writeVectorToFile(const std::vector<double> &vec,
                       const std::string &filename) {
  std::ofstream file(filename);
  if (!file.is_open()) {
    std::cerr << "Failed to open file: " << filename << std::endl;
    return;
  }

  if (vec.size() > 0) {
    file << vec[0];
    for (size_t i = 1; i < vec.size(); ++i) {
      file << ", " << vec[i];
    }
  }

  file.close();
  if (file.good()) {
    std::cout << "Vector successfully written to " << filename << std::endl;
  } else {
    std::cerr << "Error occurred while writing to the file." << std::endl;
  }
}

std::vector<double> loadVectorFromFile(const std::string &filename) {
  std::vector<double> vec;
  std::ifstream file(filename);

  if (!file.is_open()) {
    std::cerr << "Failed to open file: " << filename << std::endl;
    return vec;
  }

  std::string line;
  if (std::getline(file, line)) { // Read the single line
    std::stringstream ss(line);
    std::string value;
    while (std::getline(ss, value, ',')) { // Split by comma
      try {
        vec.push_back(std::stod(value)); // Convert to double
      } catch (const std::invalid_argument &e) {
        std::cerr << "Invalid number found: " << value << std::endl;
      }
    }
  }

  file.close();
  return vec;
}

void save_sparse_matrix(const Eigen::SparseMatrix<double> &matrix,
                        const std::string &filename) {
  // Using the Matrix Market format to save the sparse matrix
  std::ofstream file(filename);
  if (file.is_open()) {
    Eigen::SparseMatrix<double> matrix_transposed =
        matrix.transpose(); // Transpose to make saving easier
    matrix_transposed
        .makeCompressed();     // Ensure the matrix is in compressed form
    file << matrix_transposed; // Write the matrix to file
    file.close();
    std::cout << "Sparse matrix saved successfully to " << filename
              << std::endl;
  } else {
    std::cerr << "Error opening file for writing!" << std::endl;
  }
}

void save_points(const std::vector<std::array<double, 3>> &points,
                 const std::string &filename) {
  std::ofstream file(filename);

  if (file.is_open()) {
    for (const auto &point : points) {
      file << point[0] << " " << point[1] << " " << point[2] << "\n";
    }
    file.close();
    std::cout << "Points saved successfully to " << filename << std::endl;
  } else {
    std::cerr << "Error opening file for writing!" << std::endl;
  }
}

std::vector<std::array<double, 3>> load_points(const std::string &filename) {
  std::vector<std::array<double, 3>> points;
  std::ifstream file(filename);

  if (file.is_open()) {
    std::string line;
    while (std::getline(file, line)) {
      std::istringstream iss(line);
      std::array<double, 3> point;
      if (iss >> point[0] >> point[1] >> point[2]) {
        points.push_back(point);
      } else {
        std::cerr << "Invalid line in file: " << line << std::endl;
      }
    }
    file.close();
  } else {
    std::cerr << "Error opening file for reading!" << std::endl;
  }

  return points;
}
