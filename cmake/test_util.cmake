################################################################################
# Utility functions for creating tests

if(TUVX_ENABLE_MEMCHECK)
  find_program(MEMORYCHECK_COMMAND "valgrind")
endif()

################################################################################
# impose that one test runs after another so that we can safely test in parallel

function(add_test_dependency run_second run_first)
    # add dependency between two tests
    # https://stackoverflow.com/a/66931930/5217293
    set_tests_properties(${run_first}  PROPERTIES FIXTURES_SETUP    f_${run_first})
    set_tests_properties(${run_second} PROPERTIES FIXTURES_REQUIRED f_${run_first})
endfunction(add_test_dependency)

################################################################################
# build and add a standard test (one linked to the tuvx library)

function(create_standard_test)
  set(prefix TEST)
  set(singleValues NAME WORKING_DIRECTORY)
  set(multiValues SOURCES)
  include(CMakeParseArguments)
  cmake_parse_arguments(${prefix} " " "${singleValues}" "${multiValues}" ${ARGN})
  add_executable(test_${TEST_NAME} ${TEST_SOURCES})
  set_target_properties(test_${TEST_NAME} PROPERTIES LINKER_LANGUAGE Fortran)
  target_link_libraries(test_${TEST_NAME} PUBLIC musica::tuvx tuvx_test_utils GTest::gtest_main)
  if(TUVX_ENABLE_OPENMP)
    target_link_libraries(test_${TEST_NAME} PUBLIC OpenMP::OpenMP_Fortran)
  endif()
  if(NOT DEFINED TEST_WORKING_DIRECTORY)
    set(TEST_WORKING_DIRECTORY "${CMAKE_BINARY_DIR}")
  endif()
  if(TUVX_ENABLE_LAPACK)
    target_link_libraries(test_${TEST_NAME} PUBLIC ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
  endif()
  add_tuvx_test(${TEST_NAME} test_${TEST_NAME} "" ${TEST_WORKING_DIRECTORY})
endfunction(create_standard_test)

################################################################################
# build and add a standard test (one linked to the tuvx library)

function(create_standard_cxx_test)
  set(prefix TEST)
  set(optionalValues SKIP_MEMCHECK)
  set(singleValues NAME WORKING_DIRECTORY)
  set(multiValues SOURCES LIBRARIES)

  include(CMakeParseArguments)
  cmake_parse_arguments(${prefix} "${optionalValues}" "${singleValues}" "${multiValues}" ${ARGN})

  add_executable(test_${TEST_NAME} ${TEST_SOURCES})
  target_link_libraries(test_${TEST_NAME} PUBLIC musica::tuvx GTest::gtest_main)

  if(TUVX_ENABLE_LAPACK)
    target_include_directories(test_${TEST_NAME} PUBLIC ${LAPACK_INCLUDE_DIRS})
    target_link_libraries(test_${TEST_NAME} PUBLIC LAPACK::LAPACK ${LAPACKE_LIBRARIES})
  endif()

  # link additional libraries
  foreach(library ${TEST_LIBRARIES})
    target_link_libraries(test_${TEST_NAME} PUBLIC ${library})
  endforeach()

  if(NOT DEFINED TEST_WORKING_DIRECTORY)
    set(TEST_WORKING_DIRECTORY "${CMAKE_BINARY_DIR}")
  endif()

  add_tuvx_test(${TEST_NAME} test_${TEST_NAME} "" ${TEST_WORKING_DIRECTORY} ${TEST_SKIP_MEMCHECK})
endfunction(create_standard_cxx_test)

################################################################################
# Add a test

function(add_tuvx_test test_name test_binary test_args working_dir)
  if(TUVX_ENABLE_MPI)
    add_test(NAME ${test_name}
      COMMAND mpirun -v -np 2 ${CMAKE_BINARY_DIR}/${test_binary} ${test_args}
             WORKING_DIRECTORY ${working_dir})
  else()
    add_test(NAME ${test_name}
             COMMAND ${test_binary} ${test_args}
             WORKING_DIRECTORY ${working_dir})
  endif()
  set(MEMORYCHECK_COMMAND_OPTIONS "--error-exitcode=1 --trace-children=yes --leak-check=full -s --gen-suppressions=all ${MEMCHECK_SUPPRESS}")
  set(memcheck "${MEMORYCHECK_COMMAND} ${MEMORYCHECK_COMMAND_OPTIONS}")
  separate_arguments(memcheck)
  if(TUVX_ENABLE_MPI AND MEMORYCHECK_COMMAND AND TUVX_ENABLE_MEMCHECK)
    add_test(NAME memcheck_${test_name}
      COMMAND mpirun -v -np 2 ${memcheck} ${CMAKE_BINARY_DIR}/${test_binary} ${test_args}
             WORKING_DIRECTORY ${working_dir})
  elseif(MEMORYCHECK_COMMAND AND TUVX_ENABLE_MEMCHECK)
    add_test(NAME memcheck_${test_name}
             COMMAND ${memcheck} ${CMAKE_BINARY_DIR}/${test_binary} ${test_args}
             WORKING_DIRECTORY ${working_dir})
  endif()
endfunction(add_tuvx_test)

################################################################################
# Setup regression tests. Add dependencies between each regression test and its 
# memcheck test. Also add a dependence with any previous tests. Becuase TUV-x
# outputs to the same location, concurrent runs of the standalone tool that
# depend on the output must run in serial

function(add_regression_test test_name command memcheck_command)
  add_test(NAME ${test_name} COMMAND ${command} WORKING_DIRECTORY ${CMAKE_BINARY_DIR})

  if(MEMORYCHECK_COMMAND AND TUVX_ENABLE_MEMCHECK)
    add_test(NAME memcheck_${test_name} COMMAND ${memcheck_command} WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
  endif()

endfunction(add_regression_test)

################################################################################
# Link tuv-x to a test and add it to the suite as a bash script

macro(add_std_test_script test_name script_path)
  target_include_directories(${test_name} PUBLIC ${CMAKE_BINARY_DIR}/src)
  target_link_libraries(${test_name} PUBLIC musica::tuvx)
  if(TUVX_ENABLE_OPENMP)
    target_link_libraries(${test_name} PUBLIC OpenMP::OpenMP_Fortran)
  endif()
  add_test(NAME ${test_name} COMMAND ${script_path})
endmacro(add_std_test_script)

################################################################################
