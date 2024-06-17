function(create_standard_cxx_benchmark)
  set(prefix BENCHMARK)
  set(optionalValues SKIP_MEMCHECK)

  include(CMakeParseArguments)
  cmake_parse_arguments(${prefix} "${optionalValues}" "${singleValues}"
                        "${multiValues}" ${ARGN})

  add_executable(benchmark_${BENCHMARK_NAME} ${BENCHMARK_SOURCES})
  target_include_directories(benchmark_${BENCHMARK_NAME}
                             PUBLIC ${LAPACK_INCLUDE_DIRS})
  target_link_libraries(benchmark_${BENCHMARK_NAME}
                        PUBLIC LAPACK::LAPACK musica::tuvx benchmark::benchmark)

  # add_tuvx_test(${BENCHMARK_NAME} benchmark_${BENCHMARK_NAME} ""
  # ${BENCHMARK_WORKING_DIRECTORY} ${BENCHMARK_SKIP_MEMCHECK})
endfunction(create_standard_cxx_benchmark)
