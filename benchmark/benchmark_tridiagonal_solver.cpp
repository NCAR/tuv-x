
#include <tuvx/linear_algebra/linear_algebra.hpp>

#include <benchmark/benchmark.h>

#include <random>
#include <vector>

#ifdef TUVX_COMPILE_WITH_INTEL
  #include <mkl_lapacke.h>
#elif TUVX_COMPILE_WITH_GCC
  #include <lapacke.h>
#endif

const std::size_t SYSTEM_SIZE = 1e6;

const bool MAKE_DIAGONALLY_DOMINANT = true;

const unsigned RANDOM_NUMBER_SEED = 1;

/// @brief This function benchmarks the lapacke tridiagonal matrix solver for single precision
/// @param state Benchmarking argument
static void BM_LAPACKE_SINGLE_PRECISISON(benchmark::State& state)
{
  std::vector<float> x(SYSTEM_SIZE);
  std::vector<float> b(SYSTEM_SIZE);
  tuvx::TridiagonalMatrix<float> A(SYSTEM_SIZE);
  std::mt19937 random_device(RANDOM_NUMBER_SEED);
  for (auto _ : state)
  {
    state.PauseTiming();
    tuvx::FillRandom<float>(A, RANDOM_NUMBER_SEED, MAKE_DIAGONALLY_DOMINANT);
    tuvx::FillRandom<float>(x, RANDOM_NUMBER_SEED);
    b = tuvx::Dot<float>(A, x);
    state.ResumeTiming();

    LAPACKE_sgtsv(
        LAPACK_ROW_MAJOR,
        SYSTEM_SIZE,
        1,
        A.lower_diagonal_.data(),
        A.main_diagonal_.data(),
        A.upper_diagonal_.data(),
        b.data(),
        1);
  }
}

/// @brief This function benchmarks the lapacke tridiagonal matrix solver for double precision
/// @param state Benchmarking argument
static void BM_LAPACKE_DOUBLE_PRECISISON(benchmark::State& state)
{
  std::vector<double> x(SYSTEM_SIZE);
  std::vector<double> b(SYSTEM_SIZE);
  tuvx::TridiagonalMatrix<double> A(SYSTEM_SIZE);

  std::mt19937 random_device(RANDOM_NUMBER_SEED);
  // Perform setup here
  for (auto _ : state)
  {
    state.PauseTiming();
    tuvx::FillRandom<double>(A, RANDOM_NUMBER_SEED, MAKE_DIAGONALLY_DOMINANT);
    tuvx::FillRandom<double>(x, RANDOM_NUMBER_SEED);
    b = tuvx::Dot<double>(A, x);
    state.ResumeTiming();

    LAPACKE_dgtsv(
        LAPACK_ROW_MAJOR,
        SYSTEM_SIZE,
        1,
        A.lower_diagonal_.data(),
        A.main_diagonal_.data(),
        A.upper_diagonal_.data(),
        b.data(),
        1);
  }
}

/// @brief This function benchmarks the tuvx tridiagonal matrix solver for single precision
/// @param state Benchmarking argument
static void BM_TUVX_DOUBLE_PRECISISON(benchmark::State& state)
{
  std::vector<double> x(SYSTEM_SIZE);
  std::vector<double> b(SYSTEM_SIZE);
  tuvx::TridiagonalMatrix<double> A(SYSTEM_SIZE);
  std::vector<double> x_approx(SYSTEM_SIZE);

  for (auto _ : state)
  {
    state.PauseTiming();
    tuvx::FillRandom<double>(A, RANDOM_NUMBER_SEED, MAKE_DIAGONALLY_DOMINANT);
    tuvx::FillRandom<double>(x, RANDOM_NUMBER_SEED);
    state.ResumeTiming();
    tuvx::Solve<double>(A, b);
  }
}

/// @brief This function benchmarks the tuvx tridiagonal matrix solver for double precision
/// @param state Benchmarking argument
static void BM_TUVX_SINGLE_PRECISISON(benchmark::State& state)
{
  std::vector<float> x(SYSTEM_SIZE);
  std::vector<float> b(SYSTEM_SIZE);
  tuvx::TridiagonalMatrix<float> A(SYSTEM_SIZE);
  std::vector<float> x_approx(SYSTEM_SIZE);

  // Perform setup here
  for (auto _ : state)
  {
    state.PauseTiming();
    tuvx::FillRandom<float>(A, RANDOM_NUMBER_SEED, MAKE_DIAGONALLY_DOMINANT);
    tuvx::FillRandom<float>(x, RANDOM_NUMBER_SEED);
    b = tuvx::Dot<float>(A, x);
    state.ResumeTiming();
    tuvx::Solve<float>(A, b);
  }
}

/// @brief Register the functions defined above as a benchmark
BENCHMARK(BM_LAPACKE_DOUBLE_PRECISISON);
BENCHMARK(BM_LAPACKE_SINGLE_PRECISISON);
BENCHMARK(BM_TUVX_DOUBLE_PRECISISON);
BENCHMARK(BM_TUVX_SINGLE_PRECISISON);

/// @brief Run all benchmarks
BENCHMARK_MAIN();
