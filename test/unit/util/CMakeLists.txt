################################################################################
# Test utilities

include(test_util)

################################################################################
# Utility tests

create_standard_test(NAME util_array SOURCES array.F90)

create_standard_cxx_test(NAME util_array2d SOURCES array2d.cpp)

create_standard_cxx_test(NAME util_array3d SOURCES array3d.cpp)

create_standard_test(NAME util_assert SOURCES assert.F90)
add_executable(util_assert_failure assert.F90)
set_target_properties(util_assert_failure PROPERTIES LINKER_LANGUAGE Fortran)
add_std_test_script(util_assert_failure assert.sh)

create_standard_test(NAME util_config SOURCES config.F90)

create_standard_test(NAME util_map SOURCES map.F90)
add_executable(util_map_failure map.F90)
set_target_properties(util_map_failure PROPERTIES LINKER_LANGUAGE Fortran)
add_std_test_script(util_map_failure map.sh)

create_standard_test(NAME util_mpi SOURCES mpi.F90)

create_standard_test(NAME util_string SOURCES string.F90)
add_executable(util_string_failure string.F90)
set_target_properties(util_string_failure PROPERTIES LINKER_LANGUAGE Fortran)
add_std_test_script(util_string_failure string.sh)

add_subdirectory(io)

################################################################################
