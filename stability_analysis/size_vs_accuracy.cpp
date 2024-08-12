
#include <tuvx/linear_algebra/linear_algebra.hpp>

#include <benchmark/benchmark.h>

#include <array>
#include <fstream>
#include <iostream>
#include <random>
#include <vector>

#ifdef TUVX_COMPILE_WITH_INTEL
  #include <mkl_lapacke.h>
#elif TUVX_COMPILE_WITH_GCC
  #include <lapacke.h>
#endif
const double TOL_DP = std::numeric_limits<double>::epsilon();
const float TOL_SP = std::numeric_limits<float>::epsilon();

const std::size_t NUMBER_OF_RUNS = 1000;
const std::size_t SYSTEM_SIZE = 10000;
const bool MAKE_DIAGONALLY_DOMINANT = false;

const unsigned RANDOM_NUMBER_SEED = 1;

void TestLapackeNonDiagonallyDominant()
{
  std::size_t n_times_failed = 0;
  double error_run = 0;
  double error = 0;
  std::array<std::size_t, 7> sizes{ 10, 100, 1000, 10000, 100000, 1000000, 10000000 };

  std::ofstream data("data_lapacke.csv");

  for (auto s : sizes)
  {
    for (std::size_t j = 0; j < NUMBER_OF_RUNS; j++)
    {
      std::vector<double> x(s);
      std::vector<double> b(s);
      tuvx::TridiagonalMatrix<double> A(s);

      tuvx::FillRandom<double>(A, RANDOM_NUMBER_SEED, MAKE_DIAGONALLY_DOMINANT);
      tuvx::FillRandom<double>(x, RANDOM_NUMBER_SEED);

      b = tuvx::Dot<double>(A, x);

      LAPACKE_dgtsv(
          LAPACK_ROW_MAJOR, s, 1, A.lower_diagonal_.data(), A.main_diagonal_.data(), A.upper_diagonal_.data(), b.data(), 1);

      error_run = tuvx::ComputeError<double>(x, b);
      error += error_run;
    }

    error /= NUMBER_OF_RUNS;
    data << s << "," << error << "\n";
  }
  data.close();
}

void TestTuvxNonDiagonallyDominant()
{
  int n_times_failed = 0;
  double error_run = 0;
  double error = 0;
  std::array<std::size_t, 7> sizes{ 10, 100, 1000, 10000, 100000, 1000000, 10000000 };

  std::ofstream data;
  data.open("data_tuv.csv");

  for (auto s : sizes)
  {
    for (std::size_t j = 0; j < NUMBER_OF_RUNS; j++)
    {
      std::vector<double> x(s);
      std::vector<double> b(s);
      tuvx::TridiagonalMatrix<double> A(s);

      tuvx::FillRandom<double>(A, RANDOM_NUMBER_SEED, MAKE_DIAGONALLY_DOMINANT);
      tuvx::FillRandom<double>(x, RANDOM_NUMBER_SEED);

      b = tuvx::Dot<double>(A, x);

      tuvx::Solve(A, b);

      error_run = tuvx::ComputeError<double>(x, b);
      error += error_run;
    }

    error /= NUMBER_OF_RUNS;
    data << s << "," << error << "\n";
  }

  data.close();
}

int main(int argc, char *argv[])
{
  TestLapackeNonDiagonallyDominant();
  TestTuvxNonDiagonallyDominant();
  return 0;
}
