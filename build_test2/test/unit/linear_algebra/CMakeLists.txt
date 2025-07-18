# ##############################################################################
# Test utilities

include(test_util)

# ##############################################################################
# tests

create_standard_cxx_test(NAME error_function SOURCES test_error_function.cpp)

if(TUVX_ENABLE_LAPACK)
  create_standard_cxx_test(NAME tridiagonal_solver SOURCES test_tridiagonal_solver.cpp)
endif()
# ##############################################################################
