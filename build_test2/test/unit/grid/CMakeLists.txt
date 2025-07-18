################################################################################
# Test utilities

include(test_util)

################################################################################
# Photo-decomp tests

create_standard_test(NAME grid SOURCES test.F90)
create_standard_test(NAME grid_from_host SOURCES from_host.F90)

################################################################################
