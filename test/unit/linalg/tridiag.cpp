#include <cfloat>
#include <cstdlib>
#include <gtest/gtest.h>
#include <limits>
#include <tuvx/linalg/linalg.h>
#include <vector>

using namespace tuvx::linalg;

typedef trid_mat<double> trid_matd;
typedef std::vector<double> vecd;

typedef trid_mat<float> trid_matf;
typedef std::vector<float> vecf;

const int size = 10;
const double tol_dp = 1e-14;
const float tol_sp = 1e-6;

// Demonstrate some basic assertions.
TEST(TridiagSolveCPP, DoublePrecision) {
  vecd x;
  vecd b;
  trid_matd A;

  // initialize random matrix A and vector
  fill_rand_mat(A, size);
  fill_rand_vec(x, size);
  b = dot(A, x);

  // reconstruct x by tridiagonal solve
  vecd x_approx = tridiag_solve(A, b);

  // compute reconstruction error
  double error = 0.0;
  for (int i = 0; i < x_approx.size(); i++) {
    error += std::pow(std::abs(x[i] - x_approx[i]), 2);
  }
  error = std::sqrt(error) / x_approx.size();

  EXPECT_LE(error, tol_dp) << "DOUBLE EPSILON: "
                           << std::numeric_limits<double>::epsilon()
                           << std::endl;
}

TEST(TridiagSolveCPP, SinglePrecision) {
  vecf x;
  vecf b;
  trid_matf A;

  fill_rand_mat(A, size);
  fill_rand_vec(x, size);

  b = dot(A, x);

  vecf x_approx = tridiag_solve(A, b);
  float error = 0.0;
  for (int i = 0; i < x_approx.size(); i++) {
    error += std::pow(std::abs(x[i] - x_approx[i]), 2);
  }
  error = std::sqrt(error) / x_approx.size();

  EXPECT_LE(error, tol_sp) << "FLOAT EPSILON: "
                           << std::numeric_limits<float>::epsilon()
                           << std::endl;
}
