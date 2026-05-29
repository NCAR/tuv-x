#include <tuvx/linear_algebra/linear_algebra.hpp>

#include <gtest/gtest.h>

#include <cfloat>
#include <cstddef>
#include <cstdlib>
#include <limits>

const double TOL_DP = 10 * std::numeric_limits<double>::epsilon();
const float TOL_SP = 10 * std::numeric_limits<float>::epsilon();

const std::size_t NUMBER_OF_RUNS = 20;
const std::size_t SYSTEM_SIZE = 10;
const bool MAKE_DIAGONALLY_DOMINANT = true;

const unsigned RANDOM_NUMBER_SEED = 1;

/// @brief Tridiagonal Solver Test for single Precision Floats.
TEST(TridiagSolveTest, SinglePrecision)
{
  float error = 0;

  for (std::size_t j = 0; j < NUMBER_OF_RUNS; j++)
  {
    tuvx::Array1D<float> x(SYSTEM_SIZE);
    tuvx::Array1D<float> b(SYSTEM_SIZE);
    tuvx::TridiagonalMatrix<float> tridiag(SYSTEM_SIZE);

    tuvx::FillRandom<float>(tridiag, RANDOM_NUMBER_SEED, MAKE_DIAGONALLY_DOMINANT);
    tuvx::FillRandom<float>(x, RANDOM_NUMBER_SEED);
    b = tuvx::Dot<float>(tridiag, x);
    tuvx::Solve<float>(tridiag, b);

    error += tuvx::ComputeError<float>(x, b);
  }

  error /= NUMBER_OF_RUNS;

  EXPECT_LE(error, TOL_SP);
}

/// @brief Tridiagonal Solver Test for Double Precision Floats.
TEST(TridiagSolveTest, DoublePrecision)
{
  double error = 0;
  for (std::size_t j = 0; j < NUMBER_OF_RUNS; j++)
  {
    tuvx::Array1D<double> x(SYSTEM_SIZE);
    tuvx::Array1D<double> b(SYSTEM_SIZE);
    tuvx::TridiagonalMatrix<double> tridiag(SYSTEM_SIZE);

    tuvx::FillRandom<double>(tridiag, RANDOM_NUMBER_SEED, MAKE_DIAGONALLY_DOMINANT);
    tuvx::FillRandom<double>(x, RANDOM_NUMBER_SEED);
    b = tuvx::Dot<double>(tridiag, x);
    tuvx::Solve<double>(tridiag, b);

    error += tuvx::ComputeError<double>(x, b);
  }
  error /= NUMBER_OF_RUNS;
  EXPECT_LE(error, TOL_DP);
}
