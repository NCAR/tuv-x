add_library(tuvx_benchmark_utils OBJECT)

target_compile_features(tuvx_benchmark_utils PUBLIC cxx_std_11)

# link benchmark library to sources in this file
target_sources(tuvx_benchmark_utils PUBLIC benchmark_tridiagonal_solver.cpp)
target_link_libraries(tuvx_benchmark_utils PUBLIC musica::tuvx LAPACK::LAPACK
                                                  benchmark::benchmark)
target_include_directories(tuvx_benchmark_utils PUBLIC ${LAPACK_INCLUDE_DIRS})

add_executable(benchmark_tridiagonal_solver benchmark_tridiagonal_solver.cpp)
target_link_libraries(benchmark_tridiagonal_solver
                      PUBLIC LAPACK::LAPACK benchmark::benchmark musica::tuvx)
