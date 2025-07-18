################################################################################
# Test utilities

include(test_util)

################################################################################
# tests

create_standard_test(NAME spectral_weight_base SOURCES base.F90 )
create_standard_test(NAME spectral_weight_gaussian SOURCES gaussian.F90 )
create_standard_test(NAME spectral_weight_notch_filter SOURCES notch_filter.F90 )

################################################################################
