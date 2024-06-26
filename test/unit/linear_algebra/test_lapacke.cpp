#include <tuvx/linear_algebra/linear_algebra.hpp>

#include <gtest/gtest.h>

#include <cfloat>
#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <vector>

#ifdef TUVX_COMPILE_WITH_INTEL
  #include <mkl_lapacke.h>
#elif TUVX_COMPILE_WITH_GCC
  #include <lapacke.h>
#endif

using namespace tuvx;

typedef TridiagonalMatrix<double> trid_matd;
typedef std::vector<double> vecd;

typedef TridiagonalMatrix<float> trid_matf;
typedef std::vector<float> vecf;

const double tol_dp = std::numeric_limits<double>::epsilon();
const float tol_sp = std::numeric_limits<float>::epsilon();

const std::size_t number_of_runs = 100;
const std::size_t sizes[5] = { 1000, 10000, 100000, 1000000, 10000000 };

/// @test LAPACKE Tridiagonal Solver Test for single Precision Floats.
/// @brief Generate random tridiagonal matrix $A$ and vector $x$,
/// compute $b=A \cdot x$, and check if solution is reconstructed
/// accurately using L2 norm (single precision). Check for different
/// sizes to check consistency
TEST(LapackeTest, SinglePrecision)
{
  vecf errors(5, 0);

  for (std::size_t i = 0; i < 5; i++)
  {
    for (std::size_t j = 0; j < number_of_runs; j++)
    {
      vecf x(sizes[i]);
      vecf b(sizes[i]);
      trid_matf A(sizes[i]);

      FillRandom<float>(A);
      FillRandom<float>(x);
      b = Dot<float>(A, x);

      LAPACKE_sgtsv(
          LAPACK_ROW_MAJOR,
          sizes[i],
          1,
          A.lower_diagonal_.data(),
          A.main_diagonal_.data(),
          A.upper_diagonal_.data(),
          b.data(),
          1);

      // to be written to a file
      errors[i] += ComputeError<float>(x, b);
    }
    errors[i] /= number_of_runs;
    EXPECT_LE(errors[i], tol_sp);
  }
}

/// @test LAPACKE Tridiagonal Solver Test for double Precision Floats.
/// @brief Generate random tridiagonal matrix $A$ and vector $x$,
/// compute $b=A \cdot x$, and check if solution is reconstructed
/// accurately using L2 norm (double precision). Check for different
/// sizes to check consistency
TEST(LapackeTest, DoublePrecision)
{
  vecd errors(5, 0);
  vecd times(5, 0);
  for (std::size_t i = 0; i < 5; i++)
  {
    for (std::size_t j = 0; j < number_of_runs; j++)
    {
      vecd x(sizes[i]);
      vecd b(sizes[i]);
      trid_matd A(sizes[i]);

      FillRandom<double>(A);
      FillRandom<double>(x);
      b = Dot<double>(A, x);

      LAPACKE_dgtsv(
          LAPACK_ROW_MAJOR,
          sizes[i],
          1,
          A.lower_diagonal_.data(),
          A.main_diagonal_.data(),
          A.upper_diagonal_.data(),
          b.data(),
          1);

      errors[i] += ComputeError<double>(x, b);
    }
    errors[i] /= number_of_runs;
    EXPECT_LE(errors[i], tol_dp);
  }
}
