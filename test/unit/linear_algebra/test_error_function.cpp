
#include <tuvx/linear_algebra/linear_algebra.hpp>

#include <gtest/gtest.h>

#include <cfloat>
#include <cstddef>
#include <cstdlib>
#include <vector>

const std::size_t SIZE = 10;  // size of the system to test
const double TOL_DP = 1e-14;  // tolorance for double
const float TOL_SP = 1e-5;    // tolorance for single point floating precision

const unsigned RANDOM_NUMBER_SEED = 1;

/// @brief Test the correctness of the error function used for
/// testing the Linear approximation solvers with double precision data types.
TEST(ErrorFunctionTest, DoublePrecision)
{
  // same vector should return 0 error
  std::vector<double> x(SIZE);
  std::vector<double> x1(SIZE);
  tuvx::FillRandom<double>(x, RANDOM_NUMBER_SEED);
  x1 = x;
  double error = tuvx::ComputeError<double>(x, x1);
  EXPECT_EQ(error, 0.0);

  // L1 norm between x and (x[1]+0.1) should be 0.1/size;
  std::fill(x1.begin(), x1.end(), 1.0);
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
  std::vector<float> x(SIZE);
  std::vector<float> x1(SIZE);
  tuvx::FillRandom<float>(x, RANDOM_NUMBER_SEED);
  x1 = x;
  float error = tuvx::ComputeError<float>(x, x1);
  EXPECT_EQ(error, 0.0f);

  // L1 norm between x and (x[1]+0.1) should be 0.1/size;
  std::fill(x1.begin(), x1.end(), 1.0f);
  x1 = x;

  x[0] += 0.1f;
  error = tuvx::ComputeError(x, x1);
  EXPECT_LE(error - 0.1f / SIZE, TOL_SP);
}
