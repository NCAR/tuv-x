#include <tuvx/linear_algebra/linear_algebra.hpp>

#include <gtest/gtest.h>

#include <cfloat>
#include <cstddef>
#include <cstdlib>
#include <limits>

#ifdef TUVX_COMPILE_WITH_INTEL
  #include <mkl_lapacke.h>
#elif TUVX_COMPILE_WITH_GCC
  #include <lapacke.h>
#endif

const double TOL_DP = std::numeric_limits<double>::epsilon();
const float TOL_SP = std::numeric_limits<float>::epsilon();

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
    tuvx::TridiagonalMatrix<float> A(SYSTEM_SIZE);

    tuvx::FillRandom<float>(A, RANDOM_NUMBER_SEED, MAKE_DIAGONALLY_DOMINANT);
    tuvx::FillRandom<float>(x, RANDOM_NUMBER_SEED);
    b = tuvx::Dot<float>(A, x);
    tuvx::Solve<float>(A, b);

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
    tuvx::TridiagonalMatrix<double> A(SYSTEM_SIZE);

    tuvx::FillRandom<double>(A, RANDOM_NUMBER_SEED, MAKE_DIAGONALLY_DOMINANT);
    tuvx::FillRandom<double>(x, RANDOM_NUMBER_SEED);
    b = tuvx::Dot<double>(A, x);
    tuvx::Solve<double>(A, b);

    error += tuvx::ComputeError<double>(x, b);
  }
  error /= NUMBER_OF_RUNS;
  EXPECT_LE(error, TOL_DP);
}

/// @brief LAPACKE Tridiagonal Solver Test for single Precision Floats.
TEST(LapackeTest, SinglePrecision)
{
  float error = 0;

  for (std::size_t j = 0; j < NUMBER_OF_RUNS; j++)
  {
    tuvx::Array1D<float> x(SYSTEM_SIZE);
    tuvx::Array1D<float> b(SYSTEM_SIZE);
    tuvx::TridiagonalMatrix<float> A(SYSTEM_SIZE);

    tuvx::FillRandom<float>(A, RANDOM_NUMBER_SEED, MAKE_DIAGONALLY_DOMINANT);
    tuvx::FillRandom<float>(x, RANDOM_NUMBER_SEED);
    b = tuvx::Dot<float>(A, x);

    LAPACKE_sgtsv(
        LAPACK_ROW_MAJOR,
        SYSTEM_SIZE,
        1,
        A.lower_diagonal_.AsVector().data(),
        A.main_diagonal_.AsVector().data(),
        A.upper_diagonal_.AsVector().data(),
        b.AsVector().data(),
        1);

    error += tuvx::ComputeError<float>(x, b);
  }
  error /= NUMBER_OF_RUNS;
  EXPECT_LE(error, TOL_SP);
}

/// @brief LAPACKE Tridiagonal Solver Test for double Precision Floats.
TEST(LapackeTest, DoublePrecision)
{
  double error = 0;
  for (std::size_t j = 0; j < NUMBER_OF_RUNS; j++)
  {
    tuvx::Array1D<double> x(SYSTEM_SIZE);
    tuvx::Array1D<double> b(SYSTEM_SIZE);
    tuvx::TridiagonalMatrix<double> A(SYSTEM_SIZE);

    tuvx::FillRandom<double>(A, RANDOM_NUMBER_SEED, MAKE_DIAGONALLY_DOMINANT);
    tuvx::FillRandom<double>(x, RANDOM_NUMBER_SEED);
    b = tuvx::Dot<double>(A, x);

    LAPACKE_dgtsv(
        LAPACK_ROW_MAJOR,
        SYSTEM_SIZE,
        1,
        A.lower_diagonal_.AsVector().data(),
        A.main_diagonal_.AsVector().data(),
        A.upper_diagonal_.AsVector().data(),
        b.AsVector().data(),
        1);

    error += tuvx::ComputeError<double>(x, b);
  }
  error /= NUMBER_OF_RUNS;
  EXPECT_LE(error, TOL_DP);
}
