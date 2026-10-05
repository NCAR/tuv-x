################################################################################
# Test utilities

include(test_util)

################################################################################
# tests

create_standard_test(NAME radiative_transfer SOURCES radiative_transfer_core.F90
                                                     radiative_transfer_test.F90)

create_standard_test(NAME radiative_transfer_mpi SOURCES mpi.F90)

create_standard_test(NAME solver_from_host SOURCES solver_from_host.F90)

# windows doesn't run .sh files
if(NOT WIN32)
  add_executable(solver_from_host_failure solver_from_host.F90)
  set_target_properties(solver_from_host_failure PROPERTIES LINKER_LANGUAGE Fortran)
  target_link_libraries(solver_from_host_failure PUBLIC tuvx_test_utils)
  add_std_test_script(solver_from_host_failure solver_from_host.sh)
endif()

################################################################################
