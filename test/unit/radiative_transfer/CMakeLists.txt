################################################################################
# Test utilities

include(test_util)

################################################################################
# tests

create_standard_test(NAME radiative_transfer SOURCES radiative_transfer_core.F90
                                                     radiative_transfer_test.F90)

create_standard_test(NAME radiative_transfer_mpi SOURCES mpi.F90)

################################################################################
