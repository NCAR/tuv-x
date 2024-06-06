
#include <cfloat>
#include <cstddef>
#include <cstdlib>
#include <gtest/gtest.h>
#include <tuvx/linear_algebra/linear_algebra.hpp>
#include <vector>

const std::size_t size = 10; // size of the system to test
const int norm_order = 2;    // L2 norm for computing error
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
TEST(ErrorFunctionTest, ErrorSinglePrecision) {

  // same vector should return 0 error
  vecf x(size);
  vecf x1(size);
  FillRandom<float>(x);
  x1 = x;
  float error = ComputeError<float>(x, x1, norm_order);
  EXPECT_EQ(error, (float)0);

  // L1 norm between [0.1 0 0 0 ...] and [0 0 0 0 0] should be 0.1/size;
  std::fill(x.begin(), x.end(), (float)0);
  std::fill(x1.begin(), x1.end(), (float)0);
  x[0] = 0.1;
  error = ComputeError(x, x1, 1);
  EXPECT_LE(error - (float)0.1 / size, tol_sp);
}

/// @test Error function test
/// @brief Test the correctness of the error function used for
/// testing the Linear approximation solvers
TEST(ErrorFunctionTest, ErrorDoublePrecision) {

  // same vector should return 0 error
  vecd x(size);
  vecd x1(size);
  FillRandom<double>(x);
  x1 = x;
  double error = ComputeError<double>(x, x1, norm_order);
  EXPECT_EQ(error, (double)0);

  // L1 norm between [0.1 0 0 0 ...] and [0 0 0 0 0] should be 0.1/size;
  std::fill(x.begin(), x.end(), (double)0);
  std::fill(x1.begin(), x1.end(), (double)0);
  x[0] = 0.1;
  error = ComputeError(x, x1, 1);
  EXPECT_LE(error - (double)0.1 / size, tol_dp);
}
