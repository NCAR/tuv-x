#include <cstdlib>
#include <gtest/gtest.h>
#include <tuvx/linalg/linalg.h>
#include <vector>

using namespace tuvx::linalg;

typedef trid_mat<double> trid_matd;
typedef std::vector<double> vecd;

typedef trid_mat<float> trid_matf;
typedef std::vector<float> vecf;

const int size = 10;

// Demonstrate some basic assertions.
TEST(LinearAlgebraCPP, TridiagSolveDouble) {
  vecd x;
  vecd b;
  trid_matd A;

  double tol = 1e-6;
  fill_rand_mat(A, size);

  fill_rand_vec(x, size);
  b = dot(A, x);

  vecd x_approx = tridiag_solve(A, b);

  double error = 0.0;
  for (int i = 0; i < x_approx.size(); i++) {
    error += std::abs(x[i] - x_approx[i]);
  }
  error = error / x_approx.size();

  EXPECT_TRUE(error < tol);
}

TEST(LinearAlgebraCPP, TridiagSolveFloat) {
  vecf x;
  vecf b;
  trid_matf A;

  float tol = 1e-4;
  fill_rand_mat(A, size);

  fill_rand_vec(x, size);
  b = dot(A, x);

  vecf x_approx = tridiag_solve(A, b);

  float error = 0.0;
  for (int i = 0; i < x_approx.size(); i++) {
    error += std::abs(x[i] - x_approx[i]);
  }
  error = error / x_approx.size();

  EXPECT_TRUE(error <= tol);
}
