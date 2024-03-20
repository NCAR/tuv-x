################################################################################
# Test utilities

add_library(tuvx_test_utils unit/test_utils.F90)
target_link_libraries(tuvx_test_utils musica::tuvx)

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

string(REGEX MATCH "^[0-9]+" CMAKE_Fortran_COMPILER_MAJOR_VERSION ${CMAKE_Fortran_COMPILER_VERSION})

add_subdirectory(unit)
if(NOT ${CMAKE_Fortran_COMPILER_ID} MATCHES "NAG" AND NOT (${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU" AND ${CMAKE_Fortran_COMPILER_MAJOR_VERSION} VERSION_EQUAL 13))
  # oldtuv doesn't build with NAG, so bypass the regression tests
  add_subdirectory(oldtuv)
  add_subdirectory(regression)
endif()

################################################################################
# Run examples as tests

add_custom_target(make-tuv54-example-dir ALL COMMAND ${CMAKE_COMMAND}
                  -E make_directory ${CMAKE_BINARY_DIR}/example_tuv_5_4)
add_custom_target(link-tuv54-example-data ALL COMMAND ${CMAKE_COMMAND}
                  -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/../data ${CMAKE_BINARY_DIR}/example_tuv_5_4/data)
add_test(NAME TUV_5_4 COMMAND tuv-x ../examples/tuv_5_4.json
         WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/example_tuv_5_4)
add_custom_target(make-ts1-tsmlt-example-dir ALL COMMAND ${CMAKE_COMMAND}
                  -E make_directory ${CMAKE_BINARY_DIR}/example_ts1_tsmlt)
add_custom_target(link-ts1-tsmlt-example-data ALL COMMAND ${CMAKE_COMMAND}
                  -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/../data ${CMAKE_BINARY_DIR}/example_ts1_tsmlt/data)
add_test(NAME TS1_TSMLT COMMAND tuv-x ../examples/ts1_tsmlt.json
         WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/example_ts1_tsmlt)

add_custom_target(make-tuv54-yaml-example-dir ALL COMMAND ${CMAKE_COMMAND}
                  -E make_directory ${CMAKE_BINARY_DIR}/example_tuv_5_4_yaml)
add_custom_target(link-tuv54-yaml-example-data ALL COMMAND ${CMAKE_COMMAND}
                  -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/../data ${CMAKE_BINARY_DIR}/example_tuv_5_4_yaml/data)
add_test(NAME TUV_5_4_YAML COMMAND tuv-x ../examples/tuv_5_4.yml
         WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/example_tuv_5_4_yaml)
add_test(NAME TUV_5_4_COMPARE COMMAND python3 test/json_yaml_compare.py example_tuv_5_4 example_tuv_5_4_yaml
         WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
add_custom_target(make-ts1-tsmlt-yaml-example-dir ALL COMMAND ${CMAKE_COMMAND}
                  -E make_directory ${CMAKE_BINARY_DIR}/example_ts1_tsmlt_yaml)
add_custom_target(link-ts1-tsmlt-yaml-example-data ALL COMMAND ${CMAKE_COMMAND}
                  -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/../data ${CMAKE_BINARY_DIR}/example_ts1_tsmlt_yaml/data)
add_test(NAME TS1_TSMLT_YAML COMMAND tuv-x ../examples/ts1_tsmlt.yml
         WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/example_ts1_tsmlt_yaml)
add_test(NAME TS1_TSMLT_COMPARE COMMAND python3 test/json_yaml_compare.py example_ts1_tsmlt example_ts1_tsmlt_yaml
         WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
################################################################################
