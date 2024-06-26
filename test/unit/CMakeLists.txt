################################################################################
# Test utilities

include(test_util)

################################################################################
# Add subdirectories

add_subdirectory(cross_section)
add_subdirectory(grid)
add_subdirectory(profile)
add_subdirectory(quantum_yield)
add_subdirectory(radiative_transfer)
add_subdirectory(radiator)
add_subdirectory(spectral_weight)
add_subdirectory(tuv_doug)
add_subdirectory(linear_algebra)
add_subdirectory(util)

################################################################################
# TUV-x tests

create_standard_cxx_test(NAME cxx_grid SOURCES grid.cpp)
create_standard_cxx_test(NAME cxx_profile SOURCES profile.cpp)
create_standard_test(NAME grid_warehouse SOURCES grid_warehouse.F90)
create_standard_test(NAME heating_rates SOURCES heating_rates.F90)
create_standard_test(NAME la_sr_bands SOURCES la_sr_bands.F90 )
create_standard_test(NAME spherical_geometry SOURCES spherical_geometry.F90 )

################################################################################
