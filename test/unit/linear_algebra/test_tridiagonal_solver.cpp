#include <tuvx/linear_algebra/linear_algebra.hpp>

#include <gtest/gtest.h>

#include <cfloat>
#include <cstddef>
#include <cstdlib>
#include <iterator>
#include <limits>
#include <vector>

#ifdef TUVX_COMPILE_WITH_INTEL
  #include <mkl_lapacke.h>
#elif TUVX_COMPILE_WITH_GCC
  #include <lapacke.h>
#endif

const double tol_dp = std::numeric_limits<double>::epsilon();
const float tol_sp = std::numeric_limits<float>::epsilon();

const std::size_t number_of_runs = 20;
const std::size_t size = 10;
const bool diagonally_dominant = true;

const unsigned random_number_seed = 1;
/// @test Tridiagonal Solver Test for single Precision Floats.
/// @brief Generate random tridiagonal matrix $A$ and vector $x$,
/// compute $b=A \cdot x$, and check if solution is reconstructed
/// accurately using L2 norm (single precision). Check for different
/// sizes to check consistency
TEST(TridiagSolveTest, SinglePrecision)
{
  float error = 0;

  for (std::size_t j = 0; j < number_of_runs; j++)
  {
    std::vector<float> x(size);
    std::vector<float> b(size);
    tuvx::TridiagonalMatrix<float> A(size);

    tuvx::FillRandom<float>(A, random_number_seed, diagonally_dominant);
    tuvx::FillRandom<float>(x, random_number_seed);
    b = tuvx::Dot<float>(A, x);
    tuvx::Solve<float>(A, b);

    error += tuvx::ComputeError<float>(x, b);
  }

  error /= number_of_runs;

  EXPECT_LE(error, tol_sp);
}

/// @test Tridiagonal Solver Test for Double Precision Floats.
/// @brief Generate random tridiagonal matrix $A$ and vector $x$,
/// compute $b=A \cdot x$, and check if solution is reconstructed
/// accurately using L2 norm (double precision). Check for different
/// sizes to check consistency
TEST(TridiagSolveTest, DoublePrecision)
{
  double error = 0;
  for (std::size_t j = 0; j < number_of_runs; j++)
  {
    std::vector<double> x(size);
    std::vector<double> b(size);
    tuvx::TridiagonalMatrix<double> A(size);

    tuvx::FillRandom<double>(A, random_number_seed, diagonally_dominant);
    tuvx::FillRandom<double>(x, random_number_seed);
    b = tuvx::Dot<double>(A, x);
    tuvx::Solve<double>(A, b);

    error += tuvx::ComputeError<double>(x, b);
  }
  error /= number_of_runs;
  EXPECT_LE(error, tol_dp);
}

/// @test LAPACKE Tridiagonal Solver Test for single Precision Floats.
/// @brief Generate random tridiagonal matrix $A$ and vector $x$,
/// compute $b=A \cdot x$, and check if solution is reconstructed
/// accurately using L2 norm (single precision). Check for different
/// sizes to check consistency
TEST(LapackeTest, SinglePrecision)
{
  float error;

  for (std::size_t j = 0; j < number_of_runs; j++)
  {
    std::vector<float> x(size);
    std::vector<float> b(size);
    tuvx::TridiagonalMatrix<float> A(size);

    tuvx::FillRandom<float>(A, random_number_seed, diagonally_dominant);
    tuvx::FillRandom<float>(x, random_number_seed);
    b = tuvx::Dot<float>(A, x);

    LAPACKE_sgtsv(
        LAPACK_ROW_MAJOR, size, 1, A.lower_diagonal_.data(), A.main_diagonal_.data(), A.upper_diagonal_.data(), b.data(), 1);

    // to be written to a file
    error += tuvx::ComputeError<float>(x, b);
  }
  error /= number_of_runs;
  EXPECT_LE(error, tol_sp);
}

/// @test LAPACKE Tridiagonal Solver Test for double Precision Floats.
/// @brief Generate random tridiagonal matrix $A$ and vector $x$,
/// compute $b=A \cdot x$, and check if solution is reconstructed
/// accurately using L2 norm (double precision). Check for different
/// sizes to check consistency
TEST(LapackeTest, DoublePrecision)
{
  double error;
  for (std::size_t j = 0; j < number_of_runs; j++)
  {
    std::vector<double> x(size);
    std::vector<double> b(size);
    tuvx::TridiagonalMatrix<double> A(size);

    tuvx::FillRandom<double>(A, random_number_seed, diagonally_dominant);
    tuvx::FillRandom<double>(x, random_number_seed);
    b = tuvx::Dot<double>(A, x);

    LAPACKE_dgtsv(
        LAPACK_ROW_MAJOR, size, 1, A.lower_diagonal_.data(), A.main_diagonal_.data(), A.upper_diagonal_.data(), b.data(), 1);

    error += tuvx::ComputeError<double>(x, b);
  }
  error /= number_of_runs;
  EXPECT_LE(error, tol_dp);
}
