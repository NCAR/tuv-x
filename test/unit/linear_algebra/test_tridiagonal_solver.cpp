#include <cfloat>
#include <cstddef>
#include <cstdlib>
#include <ctime>
#include <gtest/gtest.h>
#include <limits>
#include <mkl_lapacke.h>
#include <tuvx/linear_algebra/linear_algebra.hpp>
#include <vector>

using namespace tuvx;

typedef TridiagonalMatrix<double> trid_matd;
typedef std::vector<double> vecd;

typedef TridiagonalMatrix<float> trid_matf;
typedef std::vector<float> vecf;

const std::size_t size = 100000000; // size of the system to test
const double tol_dp = 1e-14;        // tolorance for double
const float tol_sp = 1e-5; // tolorance for single point floating precision
const int norm_order = 1;  // L2 norm for computing error

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
  clock_t start = clock();
  vecd x_approx = Solve<double>(A, b);
  clock_t stop = clock();

  // reconstruct x by tridiagonal solve
  double error = ComputeError<double>(x, x_approx, norm_order);

  std::cout << "TUVX ERROR: " << error << std::endl;
  std::cout << "TUVX SPEED: " << (stop - start) / (double)CLOCKS_PER_SEC
            << std::endl;
  EXPECT_LE(error, tol_dp) << "DOUBLE EPSILON: "
                           << std::numeric_limits<double>::epsilon()
                           << std::endl;
}

TEST(TridiagSolveTest, CompareWithLAPACKE) {
  vecd x(size);
  vecd b(size);
  trid_matd A(size);

  FillRandom<double>(A);
  FillRandom<double>(x);

  b = Dot<double>(A, x);
  // vecd x_approx = Solve<double>(A, b);

  clock_t start = clock();
  LAPACKE_dgtsv(LAPACK_ROW_MAJOR, size, 1, A.lower_diagonal_.data(),
                A.main_diagonal_.data(), A.upper_diagonal_.data(), b.data(), 1);
  clock_t stop = clock();

  vecd lapacke_result = vecd(b.data(), b.data() + size);

  double lapacke_error = ComputeError<double>(lapacke_result, x, norm_order);
  // double tridiag_error = ComputeError<double>(x_approx, x, norm_order);

  std::cout << "LAPACKE ERROR: " << lapacke_error << std::endl;
  std::cout << "LAPACKE SPEED: " << (stop - start) / (double)CLOCKS_PER_SEC
            << std::endl;
  // std::cout << "TUVX-LINALG ERROR: " << tridiag_error << std::endl;

  EXPECT_LE(lapacke_error, tol_dp);
  // EXPECT_LE(tridiag_error, tol_dp);
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
