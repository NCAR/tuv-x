################################################################################
# Regression tests: radiative transfer

include(test_util)

create_standard_cxx_test(
    NAME    delta_eddington_regression
    SOURCES test_delta_eddington_regression.cpp
)

target_compile_definitions(test_delta_eddington_regression PRIVATE
    TUVX_REGRESSION_REFERENCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}/reference"
)

################################################################################
