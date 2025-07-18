################################################################################
# Test utilities

include(test_util)

################################################################################
# Photo-decomp tests

create_standard_test(NAME profile_warehouse SOURCES profile_warehouse.F90)
create_standard_test(NAME profiles SOURCES profile.F90)
create_standard_test(NAME profile_from_config SOURCES from_config.F90 )
create_standard_test(NAME profile_from_host SOURCES from_host.F90 )

################################################################################
