#include <tuvx/linear_algebra/linear_algebra.hpp>

#include <gtest/gtest.h>

#include <cfloat>
#include <chrono>
#include <cstddef>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <fstream>
#include <limits>
#include <vector>

#ifdef TUVX_COMPILE_WITH_INTEL
  #include <mkl_lapacke.h>
#elif TUVX_COMPILE_WITH_GCC
  #include <lapacke.h>
#endif

using namespace tuvx;
using std::chrono::duration;
using std::chrono::duration_cast;
using std::chrono::high_resolution_clock;
using std::chrono::milliseconds;

typedef TridiagonalMatrix<double> trid_matd;
typedef std::vector<double> vecd;

typedef TridiagonalMatrix<float> trid_matf;
typedef std::vector<float> vecf;

const double tol_dp = std::numeric_limits<double>::epsilon();
const float tol_sp = std::numeric_limits<float>::epsilon();

const std::size_t number_of_runs = 100;

/// @test LAPACKE Tridiagonal Solver Test for single Precision Floats.
/// @brief Generate random tridiagonal matrix $A$ and vector $x$,
/// compute $b=A \cdot x$, and check if solution is reconstructed
/// accurately using L2 norm (single precision). Check for different
/// sizes to check consistency
TEST(LapackeTest, SinglePrecision)
{
  std::size_t sizes[5] = { 1000, 10000, 100000, 1000000, 10000000 };
  vecf errors(5, 0);
  vecd times(5, 0);

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

      auto clock_start = high_resolution_clock::now();
      LAPACKE_sgtsv(
          LAPACK_ROW_MAJOR,
          sizes[i],
          1,
          A.lower_diagonal_.data(),
          A.main_diagonal_.data(),
          A.upper_diagonal_.data(),
          b.data(),
          1);
      auto clock_end = high_resolution_clock::now();
      duration<double, std::milli> elapsed_time = clock_end - clock_start;

      // to be written to a file
      errors[i] += ComputeError<float>(x, b);
      times[i] += elapsed_time.count();
    }
    errors[i] /= number_of_runs;
    times[i] /= number_of_runs;
    EXPECT_LE(errors[i], tol_sp);
  }

  // Open the file for output as text:
  std::ofstream file("../tool/numerical_analysis/tridiagonal_solver/lapacke_single_precision.dat");
  for (std::size_t i = 0; i < 5; i++)
  {
    file << errors[i] << " " << times[i] << std::endl;
  }
  file.close();
}

/// @test LAPACKE Tridiagonal Solver Test for double Precision Floats.
/// @brief Generate random tridiagonal matrix $A$ and vector $x$,
/// compute $b=A \cdot x$, and check if solution is reconstructed
/// accurately using L2 norm (double precision). Check for different
/// sizes to check consistency
TEST(LapackeTest, DoublePrecision)
{
  std::size_t sizes[5] = { 1000, 10000, 100000, 1000000, 10000000 };
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

      auto clock_start = high_resolution_clock::now();
      LAPACKE_dgtsv(
          LAPACK_ROW_MAJOR,
          sizes[i],
          1,
          A.lower_diagonal_.data(),
          A.main_diagonal_.data(),
          A.upper_diagonal_.data(),
          b.data(),
          1);

      auto clock_end = high_resolution_clock::now();
      duration<double, std::milli> elapsed_time = clock_end - clock_start;

      errors[i] += ComputeError<double>(x, b);
      times[i] += elapsed_time.count();
    }
    errors[i] /= number_of_runs;
    times[i] /= number_of_runs;
    EXPECT_LE(errors[i], tol_dp);
  }

  // Open the file for output as text:
  std::ofstream file("../tool/numerical_analysis/tridiagonal_solver/lapacke_double_precision.dat");
  for (std::size_t i = 0; i < 5; i++)
  {
    file << errors[i] << " " << times[i] << std::endl;
  }
  file.close();
}
