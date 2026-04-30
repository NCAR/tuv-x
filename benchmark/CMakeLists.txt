add_executable(benchmark_tridiagonal_solver benchmark_tridiagonal_solver.cpp)

target_link_libraries(benchmark_tridiagonal_solver
  PUBLIC
    benchmark::benchmark
    musica::tuvx
)
