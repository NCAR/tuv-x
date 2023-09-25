################################################################################
# Copy input files to build dir

add_custom_target(copy_doug_test_data ALL ${CMAKE_COMMAND} -E copy_directory
  ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_BINARY_DIR}/tuv_doug)

################################################################################
# library of cross section and quantum yield calculators based on Doug's
# original TUV code

add_library(tuv_doug)
target_sources(tuv_doug
  PRIVATE
    addpnt.f
    driver.F90
    gridck.f
    gridw.f
    gridz.f
    inter1.f
    inter2.f
    inter3.f
    inter4.f
)

add_subdirectory(JCALC)

target_link_libraries(tuv_doug PUBLIC musica::tuvx musica::musicacore)

################################################################################

# Tests refactored configurations based on Doug's data sets and Fortran code
add_executable(test_data_sets data_sets.F90 ../test_utils.F90)
target_link_libraries(test_data_sets PUBLIC tuv_doug)
if(ENABLE_OPENMP)
  target_link_libraries(test_data_sets PUBLIC OpenMP::OpenMP_Fortran)
endif()
add_tuvx_test(data_sets test_data_sets "" ${CMAKE_BINARY_DIR})

################################################################################
