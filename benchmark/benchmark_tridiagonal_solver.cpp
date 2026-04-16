#include <tuvx/linear_algebra/linear_algebra.hpp>

#include <benchmark/benchmark.h>

#include <random>

const std::size_t SYSTEM_SIZE = 1e6;

const bool MAKE_DIAGONALLY_DOMINANT = true;

const unsigned RANDOM_NUMBER_SEED = 1;

/// @brief Benchmarks the tuvx tridiagonal matrix solver for single precision
/// @param state Benchmarking argument
static void BM_TUVX_SINGLE_PRECISION(benchmark::State& state)
{
  tuvx::Array1D<float> x(SYSTEM_SIZE);
  tuvx::Array1D<float> b(SYSTEM_SIZE);
  tuvx::TridiagonalMatrix<float> A(SYSTEM_SIZE);

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

/// @brief Benchmarks the tuvx tridiagonal matrix solver for double precision
/// @param state Benchmarking argument
static void BM_TUVX_DOUBLE_PRECISION(benchmark::State& state)
{
  tuvx::Array1D<double> x(SYSTEM_SIZE);
  tuvx::Array1D<double> b(SYSTEM_SIZE);
  tuvx::TridiagonalMatrix<double> A(SYSTEM_SIZE);

  for (auto _ : state)
  {
    state.PauseTiming();
    tuvx::FillRandom<double>(A, RANDOM_NUMBER_SEED, MAKE_DIAGONALLY_DOMINANT);
    tuvx::FillRandom<double>(x, RANDOM_NUMBER_SEED);
    b = tuvx::Dot<double>(A, x);
    state.ResumeTiming();
    tuvx::Solve<double>(A, b);
  }
}

BENCHMARK(BM_TUVX_SINGLE_PRECISION);
BENCHMARK(BM_TUVX_DOUBLE_PRECISION);

BENCHMARK_MAIN();
