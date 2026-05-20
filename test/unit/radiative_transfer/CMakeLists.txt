################################################################################
# C++ radiative transfer unit tests

include(test_util)

################################################################################

create_standard_cxx_test(NAME spherical_geometry SOURCES test_spherical_geometry.cpp)
create_standard_cxx_test(NAME delta_eddington    SOURCES test_delta_eddington.cpp)

################################################################################
