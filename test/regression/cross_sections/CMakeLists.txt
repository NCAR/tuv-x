################################################################################
# Regression tests: cross-section algorithms

include(test_util)

create_regression_cxx_test(
    NAME    regression_cross_sections_analytic
    SOURCES analytic.cpp
)

create_regression_cxx_test(
    NAME    regression_cross_sections_tdependent
    SOURCES tdependent.cpp
)

################################################################################
