#include "mkl_lapacke.h"

#include <tuvx/linear_algebra/linear_algebra.hpp>

#include <benchmark/benchmark.h>

#include <mkl_lapack.h>

using namespace tuvx;
typedef TridiagonalMatrix<double> trid_matd;
typedef std::vector<double> vecd;

typedef TridiagonalMatrix<float> trid_matf;
typedef std::vector<float> vecf;

const std::size_t system_size = 1000;

static void BM_LAPACKE_SINGLE_PRECISISON(benchmark::State& state)
{
  vecf x(system_size);
  vecf b(system_size);
  trid_matf A(system_size);

  FillRandom<float>(A);
  FillRandom<float>(x);
  b = Dot<float>(A, x);

  // Perform setup here
  for (auto _ : state)
  {
    LAPACKE_sgtsv(
        LAPACK_ROW_MAJOR,
        system_size,
        1,
        A.lower_diagonal_.data(),
        A.main_diagonal_.data(),
        A.upper_diagonal_.data(),
        b.data(),
        1);
  }
}

static void BM_LAPACKE_DOUBLE_PRECISISON(benchmark::State& state)
{
  vecd x(system_size);
  vecd b(system_size);
  trid_matd A(system_size);

  FillRandom<double>(A);
  FillRandom<double>(x);
  b = Dot<double>(A, x);

  // Perform setup here
  for (auto _ : state)
  {
    LAPACKE_dgtsv(
        LAPACK_ROW_MAJOR,
        system_size,
        1,
        A.lower_diagonal_.data(),
        A.main_diagonal_.data(),
        A.upper_diagonal_.data(),
        b.data(),
        1);
  }
}

static void BM_TUVX_DOUBLE_PRECISISON(benchmark::State& state)
{
  vecd x(system_size);
  vecd b(system_size);
  trid_matd A(system_size);
  vecd x_approx(system_size);
  FillRandom<double>(A);
  FillRandom<double>(x);
  b = Dot<double>(A, x);

  // Perform setup here
  for (auto _ : state)
  {
    x_approx = Solve<double>(A, b);
  }
}

static void BM_TUVX_SINGLE_PRECISISON(benchmark::State& state)
{
  vecf x(system_size);
  vecf b(system_size);
  trid_matf A(system_size);
  vecf x_approx(system_size);
  FillRandom<float>(A);
  FillRandom<float>(x);
  b = Dot<float>(A, x);

  // Perform setup here
  for (auto _ : state)
  {
    x_approx = Solve<float>(A, b);
  }
}
// Register the function as a benchmark
BENCHMARK(BM_LAPACKE_DOUBLE_PRECISISON);
BENCHMARK(BM_LAPACKE_SINGLE_PRECISISON);
BENCHMARK(BM_TUVX_DOUBLE_PRECISISON);
BENCHMARK(BM_TUVX_SINGLE_PRECISISON);

// Run the benchmark
BENCHMARK_MAIN();
