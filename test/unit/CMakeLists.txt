################################################################################
# C++ unit tests

include(test_util)

add_subdirectory(linear_algebra)
add_subdirectory(util)

################################################################################
# TUV-x C++ tests

create_standard_cxx_test(NAME cxx_grid SOURCES grid.cpp)
create_standard_cxx_test(NAME cxx_profile SOURCES profile.cpp)

################################################################################
