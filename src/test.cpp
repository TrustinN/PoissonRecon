#include "BSpline.hpp"
#include "PPolynomialXd.hpp"
#include "utils/linalg.hpp"

static constexpr int Degree = 2;
static constexpr int DIM = 3;

int main() {
  // PPolynomialXD<Degree, DIM> ppXD({BSpline, BSpline, BSpline});
  // ScalarField<2> sf(BSpline);
  // std::array<double, 3> corner = {-1.5, -1.5, -1.5};
  // int max_iter = 10;
  // double step = 3.0 / max_iter;
  // for (int i = 0; i < max_iter; i++) {
  //   for (int j = 0; j < max_iter; j++) {
  //     for (int k = 0; k < max_iter; k++) {
  //       std::array<double, 3> dxyz = {i * step, j * step, k * step};
  //       std::array<double, 3> p = corner + dxyz;
  //       double eval1 = sf(p);
  //       double eval2 = ppXD(p);
  //       std::cout << "Error: " << std::abs((eval1 - eval2) / eval1)
  //                 << std::endl;
  //       std::cout << "Eval1: " << eval1 << std::endl;
  //       std::cout << "Eval2: " << eval2 << std::endl;
  //     }
  //   }
  // }

  // PPolynomialXD<Degree, DIM> ppXD({BSpline, BSpline, BSpline});
  // PPolynomialXD<Degree, DIM> ppXD2({BSpline, BSpline, BSpline});
  // PPolynomialXD<Degree, DIM> p3 = ppXD + ppXD2;
  // // std::cout << p3 << std::endl;
  // std::array<double, 3> corner = {-1.5, -1.5, -1.5};
  // int max_iter = 10;
  // double step = 3.0 / max_iter;
  // for (int i = 0; i < max_iter; i++) {
  //   for (int j = 0; j < max_iter; j++) {
  //     for (int k = 0; k < max_iter; k++) {
  //       std::array<double, 3> dxyz = {i * step, j * step, k * step};
  //       std::array<double, 3> p = corner + dxyz;
  //       double eval1 = p3(p);
  //       double eval2 = 2 * ppXD(p);
  //       std::cout << "Error: " << std::abs((eval1 - eval2) / eval1)
  //                 << std::endl;
  //       std::cout << "Eval1: " << eval1 << std::endl;
  //       std::cout << "Eval2: " << eval2 << std::endl;
  //     }
  //   }
  // }

  std::array<double, 3> corner = {-1.5, -1.5, -1.5};

  std::array<double, 3> c1 = {0, 0, 0};
  double w1 = 1.0;
  ScalarField<2> sf1(BSpline, c1, 1 / w1);

  std::array<double, 3> c2 = {0.5, -0.5, 0.5};
  double w2 = .5;
  ScalarField<2> sf2(BSpline, c2, 1 / w2);

  PPolynomialXD<Degree, DIM> pp1(sf1);
  PPolynomialXD<Degree, DIM> pp2(sf2);

  PPolynomialXD<Degree, DIM> pc = pp1 + pp2;

  int max_iter = 10;
  double step = 3.0 / max_iter;
  for (int i = 0; i < max_iter; i++) {
    for (int j = 0; j < max_iter; j++) {
      for (int k = 0; k < max_iter; k++) {
        std::array<double, 3> dxyz = {i * step, j * step, k * step};
        std::array<double, 3> p = corner + dxyz;
        double eval1 = pc(p);
        double eval2 = pp1(p) + pp2(p);
        std::cout << "Error: " << std::abs((eval1 - eval2) / eval1)
                  << std::endl;
        std::cout << "Eval1: " << eval1 << std::endl;
        std::cout << "Eval2: " << eval2 << std::endl;
      }
    }
  }

  return 0;
}
