add_executable(size_vs_accuracy size_vs_accuracy.cpp)
target_include_directories(
  size_vs_accuracy PUBLIC ${OpenBLAS_INCLUDE_DIRS}
                                      ${LAPACK_INCLUDE_DIRS})

target_link_libraries(size_vs_accuracy
  PUBLIC 
    LAPACK::LAPACK 
    ${LAPACKE_LIBRARIES}
    benchmark::benchmark 
    musica::tuvx
)
