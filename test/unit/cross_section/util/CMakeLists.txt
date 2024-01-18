################################################################################
# Test utilities

include(test_util)

################################################################################
# Cross section utility tests

create_standard_test( NAME temperature_parameterization_taylor_series
                      SOURCES temperature_parameterization_taylor_series.F90 )

################################################################################