################################################################################
# Test utilities

add_library(tuvx_test_utils unit/test_utils.F90)
target_link_libraries(tuvx_test_utils musica::musicacore)

set_target_properties(tuvx_test_utils PROPERTIES
  Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/test_include
)

target_include_directories(tuvx_test_utils 
  PUBLIC  
    ${PROJECT_BINARY_DIR}/test_include
)

include(test_util)

################################################################################
# Copy test data

add_custom_target(copy_test_data ALL ${CMAKE_COMMAND} -E copy_directory
  ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_BINARY_DIR}/test)

################################################################################
# Add subdirectories

add_subdirectory(unit)
if(NOT ${CMAKE_Fortran_COMPILER_ID} MATCHES "NAG")
  # oldtuv doesn't build with NAG, so bypass the regression tests
  add_subdirectory(oldtuv)
  add_subdirectory(regression)
endif()

################################################################################
# Run examples as tests

add_test(NAME full_example COMMAND tuv-x examples/full_config.json
         WORKING_DIRECTORY ${CMAKE_BINARY_DIR})

################################################################################
