add_executable(benchmark_tridiagonal_solver benchmark_tridiagonal_solver.cpp)
target_include_directories(
  benchmark_tridiagonal_solver PUBLIC ${OpenBLAS_INCLUDE_DIRS}
                                      ${LAPACK_INCLUDE_DIRS})

target_link_libraries(benchmark_tridiagonal_solver
                      PUBLIC LAPACK::LAPACK benchmark::benchmark musica::tuvx)
