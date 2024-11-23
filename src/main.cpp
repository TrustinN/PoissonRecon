#include "basis.hpp"
#include "utils/io.hpp"
#include <iostream>

std::ostream &operator<<(std::ostream &os, const std::array<double, 27> &arr) {
  std::string indent = std::string(3, ' ');
  os << "[" << std::endl << indent << arr[0];

  for (int i = 1; i < 27; i++) {
    os << ", ";
    if (i % 3 == 0) {
      os << std::endl;
    }
    os << indent << arr[i];
  };
  return os << std::endl << "]";
};

int main() {
  intVField v_field;
  std::cout << v_field.int_field_x << std::endl;
  std::cout << v_field.int_field_y << std::endl;
  std::cout << v_field.int_field_z << std::endl;
  return 0;
}
