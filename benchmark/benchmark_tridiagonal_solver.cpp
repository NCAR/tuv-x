
#include <tuvx/linear_algebra/linear_algebra.hpp>

#include <benchmark/benchmark.h>

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

const std::size_t system_size = 1e6;

const bool diagonally_dominant = false;

static void BM_LAPACKE_SINGLE_PRECISISON(benchmark::State& state)
{
  vecf x(system_size);
  vecf b(system_size);
  trid_matf A(system_size);

  for (auto _ : state)
  {
    state.PauseTiming();
    FillRandom<float>(A, diagonally_dominant);
    FillRandom<float>(x);
    b = Dot<float>(A, x);
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
  vecd x(system_size);
  vecd b(system_size);
  trid_matd A(system_size);

  // Perform setup here
  for (auto _ : state)
  {
    state.PauseTiming();
    FillRandom<double>(A, diagonally_dominant);
    FillRandom<double>(x);
    b = Dot<double>(A, x);
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
  vecd x(system_size);
  vecd b(system_size);
  trid_matd A(system_size);
  vecd x_approx(system_size);

  for (auto _ : state)
  {
    state.PauseTiming();
    FillRandom<double>(A, diagonally_dominant);
    FillRandom<double>(x);
    state.ResumeTiming();
    Solve<double>(A, b);
  }
}

static void BM_TUVX_SINGLE_PRECISISON(benchmark::State& state)
{
  vecf x(system_size);
  vecf b(system_size);
  trid_matf A(system_size);
  vecf x_approx(system_size);

  // Perform setup here
  for (auto _ : state)
  {
    state.PauseTiming();
    FillRandom<float>(A, diagonally_dominant);
    FillRandom<float>(x);
    b = Dot<float>(A, x);
    state.ResumeTiming();
    Solve<float>(A, b);
  }
}

// Register the function as a benchmark
BENCHMARK(BM_LAPACKE_DOUBLE_PRECISISON);
BENCHMARK(BM_LAPACKE_SINGLE_PRECISISON);
BENCHMARK(BM_TUVX_DOUBLE_PRECISISON);
BENCHMARK(BM_TUVX_SINGLE_PRECISISON);

// Run the benchmark
BENCHMARK_MAIN();
