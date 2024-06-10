
#include <cfloat>
#include <cstddef>
#include <cstdlib>
#include <gtest/gtest.h>
#include <tuvx/linear_algebra/linear_algebra.hpp>
#include <vector>

const std::size_t size = 10; // size of the system to test
const double tol_dp = 1e-14; // tolorance for double
const float tol_sp = 1e-5;   // tolorance for single point floating precision

using namespace tuvx;

typedef TridiagonalMatrix<double> trid_matd;
typedef std::vector<double> vecd;

typedef TridiagonalMatrix<float> trid_matf;
typedef std::vector<float> vecf;

/// @test Error function test
/// @brief Test the correctness of the error function used for
/// testing the Linear approximation solvers
TEST(ErrorFunctionTest, DoublePrecision) {

  // same vector should return 0 error
  vecd x(size);
  vecd x1(size);
  FillRandom<double>(x);
  x1 = x;
  double error = ComputeError<double>(x, x1);
  EXPECT_EQ(error, 0.0);

  // L1 norm between x and (x[1]+0.1) should be 0.1/size;
  std::fill(x1.begin(), x1.end(), 1.0);
  x1 = x;

  x[0] += 0.1;
  error = ComputeError(x, x1);
  EXPECT_LE(error - 0.1 / size, tol_sp);
}

/// @test Error function test
/// @brief Test the correctness of the error function used for
/// testing the Linear approximation solvers
TEST(ErrorFunctionTest, SinglePrecision) {
  // same vector should return 0 error
  vecf x(size);
  vecf x1(size);
  FillRandom<float>(x);
  x1 = x;
  float error = ComputeError<float>(x, x1);
  EXPECT_EQ(error, 0.0f);

  // L1 norm between x and (x[1]+0.1) should be 0.1/size;
  std::fill(x1.begin(), x1.end(), 1.0f);
  x1 = x;

  x[0] += 0.1f;
  error = ComputeError(x, x1);
  EXPECT_LE(error - 0.1f / size, tol_sp);
}
