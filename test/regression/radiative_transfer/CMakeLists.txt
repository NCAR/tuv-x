################################################################################
# Regression tests: radiative transfer

include(test_util)

create_standard_cxx_test(
    NAME    regression_delta_eddington
    SOURCES delta_eddington.cpp
)

target_compile_definitions(test_regression_delta_eddington PRIVATE
    TUVX_REGRESSION_REFERENCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}/reference"
)

create_standard_cxx_test(
    NAME    regression_spherical_geometry
    SOURCES spherical_geometry.cpp
)

target_compile_definitions(test_regression_spherical_geometry PRIVATE
    TUVX_REGRESSION_REFERENCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}/reference"
)

################################################################################
