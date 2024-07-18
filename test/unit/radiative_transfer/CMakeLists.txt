################################################################################
# Test utilities

include(test_util)

################################################################################
# tests

create_standard_test(NAME radiative_transfer SOURCES radiative_transfer_core.F90
                                                     rradiative_transfer_test.F90
                                                     test_source_functions.cpp)

create_standard_test(NAME radiative_transfer_mpi SOURCES mpi.F90)

################################################################################
