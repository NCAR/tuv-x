################################################################################
# C++ unit tests

include(test_util)

add_subdirectory(linear_algebra)
add_subdirectory(util)
add_subdirectory(radiative_transfer)
add_subdirectory(transforms)
add_subdirectory(interpolate)

################################################################################
# TUV-x tests

create_standard_cxx_test(NAME cxx_grid SOURCES grid.cpp)
create_standard_cxx_test(NAME cxx_profile SOURCES profile.cpp)

################################################################################
