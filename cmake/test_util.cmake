################################################################################
# Utility functions for creating C++ tests

if(TUVX_ENABLE_MEMCHECK)
  find_program(MEMORYCHECK_COMMAND "valgrind")
endif()

################################################################################
# Impose that one test runs after another so that we can safely test in parallel

function(add_test_dependency run_second run_first)
    set_tests_properties(${run_first}  PROPERTIES FIXTURES_SETUP    f_${run_first})
    set_tests_properties(${run_second} PROPERTIES FIXTURES_REQUIRED f_${run_first})
endfunction(add_test_dependency)

################################################################################
# Build and add a C++ test

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

  foreach(library ${TEST_LIBRARIES})
    target_link_libraries(test_${TEST_NAME} PUBLIC ${library})
  endforeach()

  if(NOT DEFINED TEST_WORKING_DIRECTORY)
    set(TEST_WORKING_DIRECTORY "${CMAKE_BINARY_DIR}")
  endif()

  add_tuvx_test(${TEST_NAME} test_${TEST_NAME} "" ${TEST_WORKING_DIRECTORY} ${TEST_SKIP_MEMCHECK})
endfunction(create_standard_cxx_test)

################################################################################
# Add a test (with optional memcheck)

function(add_tuvx_test test_name test_binary test_args working_dir)
  add_test(NAME ${test_name}
           COMMAND ${test_binary} ${test_args}
           WORKING_DIRECTORY ${working_dir})

  set(MEMORYCHECK_COMMAND_OPTIONS "--error-exitcode=1 --trace-children=yes --leak-check=full -s --gen-suppressions=all ${MEMCHECK_SUPPRESS}")
  set(memcheck "${MEMORYCHECK_COMMAND} ${MEMORYCHECK_COMMAND_OPTIONS}")
  separate_arguments(memcheck)

  if(MEMORYCHECK_COMMAND AND TUVX_ENABLE_MEMCHECK)
    add_test(NAME memcheck_${test_name}
             COMMAND ${memcheck} ${CMAKE_BINARY_DIR}/${test_binary} ${test_args}
             WORKING_DIRECTORY ${working_dir})
  endif()
endfunction(add_tuvx_test)

################################################################################
