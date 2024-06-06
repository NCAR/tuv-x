#include <cfloat>
#include <cstddef>
#include <cstdlib>
#include <gtest/gtest.h>
#include <limits>
#include <tuvx/linear_algebra/linear_algebra.hpp>
#include <vector>

using namespace tuvx;

typedef TridiagonalMatrix<double> trid_matd;
typedef std::vector<double> vecd;

typedef TridiagonalMatrix<float> trid_matf;
typedef std::vector<float> vecf;

const std::size_t size = 10; // size of the system to test
const double tol_dp = 1e-14; // tolorance for double
const float tol_sp = 1e-6;   // tolorance for single point floating precision
const int norm_order = 2;    // L2 norm for computing error

/// @test Tridiagonal Solver Test for Double Precision Floats.
/// @brief Generate random tridiagonal matrix $A$ and vector $x$,
/// compute $b=A \cdot x$, and check if solution is reconstructed
/// accurately using L2 norm (double precision).
TEST(TridiagSolveTest, DoublePrecision) {
  vecd x(size);
  vecd b(size);
  trid_matd A(size);

  // initialize random matrix A and vector
  FillRandom<double>(A);
  FillRandom<double>(x);
  b = Dot<double>(A, x);

  // reconstruct x by tridiagonal solve
  vecd x_approx = Solve<double>(A, b);

  // reconstruct x by tridiagonal solve
  double error = ComputeError<double>(x, x_approx, norm_order);

  EXPECT_LE(error, tol_dp) << "DOUBLE EPSILON: "
                           << std::numeric_limits<double>::epsilon()
                           << std::endl;
}

/// @test Tridiagonal Solver Test for single Precision Floats.
/// @brief Generate random tridiagonal matrix $A$ and vector $x$,
/// compute $b=A \cdot x$, and check if solution is reconstructed
/// accurately using L2 norm (single precision).
TEST(TridiagSolveTest, SinglePrecision) {
  vecf x(size);
  vecf b(size);
  trid_matf A(size);

  FillRandom<float>(A);
  FillRandom<float>(x);

  b = Dot<float>(A, x);
  vecf x_approx = Solve<float>(A, b);

  float error = ComputeError<float>(x, x_approx, norm_order);

  EXPECT_LE(error, tol_sp) << "FLOAT EPSILON: "
                           << std::numeric_limits<float>::epsilon()
                           << std::endl;
}
