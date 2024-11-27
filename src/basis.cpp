#include "basis.hpp"
#include <cassert>

constexpr static double MAX_OUTER_T = 1.5;
constexpr static double MAX_INNER_T = 0.5;
constexpr static double MIN_INNER_T = -0.5;
constexpr static double MIN_OUTER_T = -1.5;
constexpr static double EPSILON = 1e-10;

int sign(double x) { return (x > 0) - (x < 0); }

// t is centered at f1
double integral_f1_f2(int t) {
  switch (t) {
  case 0:
    return 0.00833333;
  case 1:
    return 0.55;
  case 2:
    return 0.00833333;
  default:
    assert(false);
  }
}

// t is centered at f1
double integral_df1_f2(int t) {
  switch (t) {
  case 0:
    return -0.0416667;
  case 1:
    return 0;
  case 2:
    return 0.0416667;
  default:
    assert(false);
  }
}

// t is centered at f1
double integral_d2f1_f2(int t) {
  switch (t) {
  case 0:
    return 0.166667;
  case 1:
    return -1;
  case 2:
    return 0.166667;
  default:
    assert(false);
  }
}

divergenceField::divergenceField() {
  for (int k = 0; k < 3; k++) {
    for (int j = 0; j < 3; j++) {
      for (int i = 0; i < 3; i++) {
        double entry =
            integral_df1_f2(i) * integral_f1_f2(j) * integral_f1_f2(k);

        int idx1 = i + 3 * j + 9 * k;
        x_field[idx1] = entry;

        int idx2 = k + 3 * i + 9 * j;
        y_field[idx2] = entry;

        int idx3 = j + 3 * k + 9 * i;
        z_field[idx3] = entry;
      }
    }
  }
};

laplacianField::laplacianField() {
  for (int k = 0; k < 3; k++) {
    for (int j = 0; j < 3; j++) {
      for (int i = 0; i < 3; i++) {
        double entry =
            integral_d2f1_f2(i) * integral_f1_f2(j) * integral_f1_f2(k);

        int idx1 = i + 3 * j + 9 * k;
        x_field[idx1] = entry;

        int idx2 = k + 3 * i + 9 * j;
        y_field[idx2] = entry;

        int idx3 = j + 3 * k + 9 * i;
        z_field[idx3] = entry;
      }
    }
  }
};
