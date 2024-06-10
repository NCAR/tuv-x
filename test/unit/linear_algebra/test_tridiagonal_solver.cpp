#include <cfloat>
#include <cstddef>
#include <cstdlib>
#include <ctime>
#include <gtest/gtest.h>
#include <limits>
#include <tuvx/linear_algebra/linear_algebra.hpp>
#include <vector>

using namespace tuvx;

typedef TridiagonalMatrix<double> trid_matd;
typedef std::vector<double> vecd;

typedef TridiagonalMatrix<float> trid_matf;
typedef std::vector<float> vecf;

const double tol_dp = std::numeric_limits<double>::epsilon();
const float tol_sp = std::numeric_limits<float>::epsilon();
const int norm_order = 2; // L2 norm for computing error

/// @test Tridiagonal Solver Test for single Precision Floats.
/// @brief Generate random tridiagonal matrix $A$ and vector $x$,
/// compute $b=A \cdot x$, and check if solution is reconstructed
/// accurately using L2 norm (single precision). Check for different
/// sizes to check consistency
TEST(TridiagSolveTest, SinglePrecision) {
  std::size_t sizes[5] = {5, 10, 1000, 100000, 10000000};
  float error = 0;
  for (std::size_t i = 0; i < 5; i++) {
    vecf x(sizes[i]);
    vecf b(sizes[i]);
    trid_matf A(sizes[i]);

    FillRandom<float>(A);
    FillRandom<float>(x);
    b = Dot<float>(A, x);

    vecf x_approx = Solve<float>(A, b);
    error = ComputeError<float>(x, x_approx, norm_order);
    EXPECT_LE(error, tol_sp);
  }
}

/// @test Tridiagonal Solver Test for Double Precision Floats.
/// @brief Generate random tridiagonal matrix $A$ and vector $x$,
/// compute $b=A \cdot x$, and check if solution is reconstructed
/// accurately using L2 norm (double precision). Check for different
/// sizes to check consistency
TEST(TridiagSolveTest, DoublePrecision) {
  std::size_t sizes[5] = {5, 10, 1000, 100000, 10000000};
  double error = 0;
  for (std::size_t i = 0; i < 5; i++) {
    vecd x(sizes[i]);
    vecd b(sizes[i]);
    trid_matd A(sizes[i]);

    FillRandom<double>(A);
    FillRandom<double>(x);
    b = Dot<double>(A, x);

    vecd x_approx = Solve<double>(A, b);
    error = ComputeError<double>(x, x_approx, norm_order);
    EXPECT_LE(error, tol_sp);
  }
}
