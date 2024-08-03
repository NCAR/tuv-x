
#include <tuvx/linear_algebra/linear_algebra.hpp>

#include <benchmark/benchmark.h>

#include <random>
#include <vector>

#ifdef TUVX_COMPILE_WITH_INTEL
  #include <mkl_lapacke.h>
#elif TUVX_COMPILE_WITH_GCC
  #include <lapacke.h>
#endif

const bool MAKE_DIAGONALLY_DOMINANT = true;

const unsigned RANDOM_NUMBER_SEED = 1;

void TestLapackeDiagonallyDominant()
{
}
