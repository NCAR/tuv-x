################################################################################
# Regression tests: spectral-weight algorithms

include(test_util)

create_regression_cxx_test(
    NAME    regression_spectral_weights
    SOURCES spectral_weights.cpp
)

################################################################################
