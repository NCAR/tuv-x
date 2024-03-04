################################################################################
# Test utilities

include(test_util)

################################################################################
# IO Utility tests

# There is a small memory leak in the NetCDF library code when
# creating a new file, so skip the memory check
# Also, since this creates a file, don't run with multiple cores
add_executable(test_util_io_netcdf netcdf.F90)
target_link_libraries(test_util_io_netcdf PUBLIC musica::tuvx)
target_include_directories(test_util_io_netcdf PUBLIC ${CMAKE_BINARY_DIR}/src)
add_test(NAME util_io_netcdf COMMAND ${CMAKE_BINARY_DIR}/test_util_io_netcdf)

################################################################################
