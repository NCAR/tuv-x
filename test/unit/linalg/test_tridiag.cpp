#include <gtest/gtest.h>
#include <tuvx/linalg/linalg.h>
#include <vector>

using namespace tuvx::linalg;
typedef trid_mat<double> trid_matd;
typedef std::vector<double> vecd;

// Demonstrate some basic assertions.
TEST(HelloTest, BasicAssertions) {
  vecd x;
  vecd b;
  trid_matd A;

  int size = 5;
  double tol = 1e-7;
  fill_rand_mat(A, size);
  fill_rand_vec(x, size);
  b = dot(A, x);

  tuvx::linalg::print_vec(x);
  tuvx::linalg::print_trid_mat(A);

  vecd x_approx = tridiag_solve(A, b);

  double error = 0.0;
  for (int i = 0; i < x_approx.size(); i++) {
    error += std::pow(x[i] - x_approx[i], 2);
  }
  error = std::sqrt(error);

  tuvx::linalg::print_vec(x_approx);
  EXPECT_TRUE(error < tol);
}
