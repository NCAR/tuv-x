################################################################################
# Test utilities

include(test_util)

################################################################################
# Phototolysis rate regression tests

create_standard_test(NAME cpp_delta_eddington
                     SOURCES delta_eddington.F90 delta_eddington.cpp)

################################################################################