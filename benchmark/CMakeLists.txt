add_executable(benchmark_tridiagonal_solver benchmark_tridiagonal_solver.cpp)

find_package(CUDA)

target_include_directories(
  benchmark_tridiagonal_solver PUBLIC ${OpenBLAS_INCLUDE_DIRS}
                                      ${LAPACK_INCLUDE_DIRS}
                                      ${CUDA_INCLUDE_DIRS}
                                      "/glade/derecho/scratch/sunjian/magma-2.8.0/include/")

target_link_libraries(benchmark_tridiagonal_solver
  PUBLIC 
    LAPACK::LAPACK 
    ${LAPACKE_LIBRARIES}
    ${CUDA_LIBRARIES}
    benchmark::benchmark 
    "/glade/derecho/scratch/sunjian/magma-2.8.0/lib/libmagma.so"
    musica::tuvx
)
