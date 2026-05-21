################################################################################
# Regression tests: radiative transfer

include(test_util)

create_regression_cxx_test(
    NAME    regression_delta_eddington
    SOURCES delta_eddington.cpp
)

create_regression_cxx_test(
    NAME    regression_spherical_geometry
    SOURCES spherical_geometry.cpp
)

################################################################################
