
#include <tuvx/linear_algebra/linear_algebra.hpp>

#include <gtest/gtest.h>

#include <cfloat>
#include <cstddef>
#include <cstdlib>

const std::size_t SIZE = 10;  // size of the system to test
const double TOL_DP = 1e-14;  // tolerance for double
const float TOL_SP = 1e-5;    // tolerance for single point floating precision

const unsigned RANDOM_NUMBER_SEED = 1;

/// @brief Test the correctness of the error function used for
/// testing the Linear approximation solvers with double precision data types.
TEST(ErrorFunctionTest, DoublePrecision)
{
  // same vector should return 0 error
  tuvx::Array1D<double> x(SIZE);
  tuvx::Array1D<double> x1(SIZE);
  tuvx::FillRandom<double>(x, RANDOM_NUMBER_SEED);
  x1 = x;
  double error = tuvx::ComputeError<double>(x, x1);
  EXPECT_EQ(error, 0.0);

  // L1 norm between x and (x[1]+0.1) should be 0.1/size;
  x1 = x;

  x[0] += 0.1;
  error = tuvx::ComputeError<double>(x, x1);
  EXPECT_LE(error - 0.1 / SIZE, TOL_SP);
}

/// @brief Test the correctness of the error function used for
/// testing the Linear approximation solvers with single precision data types.
TEST(ErrorFunctionTest, SinglePrecision)
{
  // same vector should return 0 error
  tuvx::Array1D<float> x(SIZE);
  tuvx::Array1D<float> x1(SIZE);
  tuvx::FillRandom<float>(x, RANDOM_NUMBER_SEED);
  x1 = x;
  float error = tuvx::ComputeError<float>(x, x1);
  EXPECT_EQ(error, 0.0f);

  // L1 norm between x and (x[1]+0.1) should be 0.1/size;
  x1 = x;

  x[0] += 0.1f;
  error = tuvx::ComputeError(x, x1);
  EXPECT_LE(error - 0.1f / SIZE, TOL_SP);
}
