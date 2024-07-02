
#include <tuvx/linear_algebra/linear_algebra.hpp>

#include <benchmark/benchmark.h>

#include <random>
#include <vector>

#ifdef TUVX_COMPILE_WITH_INTEL
  #include <mkl_lapacke.h>
#elif TUVX_COMPILE_WITH_GCC
  #include <lapacke.h>
#endif

const std::size_t system_size = 1e6;

const bool diagonally_dominant = true;

const unsigned random_number_seed = 1;

static void BM_LAPACKE_SINGLE_PRECISISON(benchmark::State& state)
{
  std::vector<float> x(system_size);
  std::vector<float> b(system_size);
  tuvx::TridiagonalMatrix<float> A(system_size);
  std::mt19937 random_device(random_number_seed);
  for (auto _ : state)
  {
    state.PauseTiming();
    tuvx::FillRandom<float>(A, random_number_seed, diagonally_dominant);
    tuvx::FillRandom<float>(x, random_number_seed);
    b = tuvx::Dot<float>(A, x);
    state.ResumeTiming();

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
  std::vector<double> x(system_size);
  std::vector<double> b(system_size);
  tuvx::TridiagonalMatrix<double> A(system_size);

  std::mt19937 random_device(random_number_seed);
  // Perform setup here
  for (auto _ : state)
  {
    state.PauseTiming();
    tuvx::FillRandom<double>(A, random_number_seed, diagonally_dominant);
    tuvx::FillRandom<double>(x, random_number_seed);
    b = tuvx::Dot<double>(A, x);
    state.ResumeTiming();

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
  std::vector<double> x(system_size);
  std::vector<double> b(system_size);
  tuvx::TridiagonalMatrix<double> A(system_size);
  std::vector<double> x_approx(system_size);

  for (auto _ : state)
  {
    state.PauseTiming();
    tuvx::FillRandom<double>(A, random_number_seed, diagonally_dominant);
    tuvx::FillRandom<double>(x, random_number_seed);
    state.ResumeTiming();
    tuvx::Solve<double>(A, b);
  }
}

static void BM_TUVX_SINGLE_PRECISISON(benchmark::State& state)
{
  std::vector<float> x(system_size);
  std::vector<float> b(system_size);
  tuvx::TridiagonalMatrix<float> A(system_size);
  std::vector<float> x_approx(system_size);

  // Perform setup here
  for (auto _ : state)
  {
    state.PauseTiming();
    tuvx::FillRandom<float>(A, random_number_seed, diagonally_dominant);
    tuvx::FillRandom<float>(x, random_number_seed);
    b = tuvx::Dot<float>(A, x);
    state.ResumeTiming();
    tuvx::Solve<float>(A, b);
  }
}

// Register the function as a benchmark
BENCHMARK(BM_LAPACKE_DOUBLE_PRECISISON);
BENCHMARK(BM_LAPACKE_SINGLE_PRECISISON);
BENCHMARK(BM_TUVX_DOUBLE_PRECISISON);
BENCHMARK(BM_TUVX_SINGLE_PRECISISON);

// Run the benchmark
BENCHMARK_MAIN();
