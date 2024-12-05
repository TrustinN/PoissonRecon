#include "basis.hpp"
#include "BSpline.hpp"
#include <cassert>

constexpr static double EPSILON = 1e-10;

int sign(double x) { return (x > 0) - (x < 0); }

// -------------------------------------------------------------------------------------------------//
// Fields
// -------------------------------------------------------------------------------------------------//

divergenceField::divergenceField() {
  auto bS = ScalarField<2>(BSpline);
  for (int k = 0; k < 3; k++) {
    for (int j = 0; j < 3; j++) {
      for (int i = 0; i < 3; i++) {
        std::array<double, 3> offset = {(double)i - 1, (double)j - 1,
                                        (double)k - 1};
        auto bS_div = ScalarField<2>(BSpline, offset);

        int idx = i + 3 * j + 9 * k;
        x_field[idx] = bS_div.partialDerivative(0).innerProduct(bS);

        y_field[idx] = bS_div.partialDerivative(1).innerProduct(bS);

        z_field[idx] = bS_div.partialDerivative(2).innerProduct(bS);
      }
    }
  }
};

laplacianField::laplacianField() {
  auto bS = ScalarField<2>(BSpline);
  for (int k = 0; k < 3; k++) {
    for (int j = 0; j < 3; j++) {
      for (int i = 0; i < 3; i++) {
        std::array<double, 3> offset = {(double)i - 1, (double)j - 1,
                                        (double)k - 1};
        auto bS_div = ScalarField<2>(BSpline, offset);

        int idx = i + 3 * j + 9 * k;
        x_field[idx] =
            bS_div.partialDerivative(0).partialDerivative(0).innerProduct(bS);

        y_field[idx] =
            bS_div.partialDerivative(1).partialDerivative(1).innerProduct(bS);

        z_field[idx] =
            bS_div.partialDerivative(2).partialDerivative(2).innerProduct(bS);
      }
    }
  }
};
