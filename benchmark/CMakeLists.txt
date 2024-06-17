# ##############################################################################
# Benchmark utilities

include(benchmark_utils)

# ##############################################################################
# tests

create_standard_cxx_benchmark(NAME tridiagonal_solver SOURCES
                              benchmark_tridiagonal_solver.cpp)

# ##############################################################################
