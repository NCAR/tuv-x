################################################################################
# Test utilities

include(test_util)

################################################################################
# tests

create_standard_test(NAME radiator SOURCES radiator_core.F90
                                           radiator_test.F90)
create_standard_test(NAME radiator_from_host SOURCES from_host.F90)

################################################################################
