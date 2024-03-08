################################################################################
# Photo-decomp tool source

# object library
add_library(tuvx_object OBJECT)

target_compile_features(tuvx_object INTERFACE cxx_std_11)

set_target_properties(tuvx_object PROPERTIES
  Fortran_MODULE_DIRECTORY ${TUVX_MOD_DIR}
)

target_include_directories(tuvx_object
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/../include>
    $<INSTALL_INTERFACE:${INSTALL_INCLUDE_DIR}>
)

message(INFO "lapack libraries: ${LAPACK_LIBRARIES}")

target_link_libraries(tuvx_object 
  PUBLIC
    PkgConfig::netcdff
    yaml-cpp::yaml-cpp
    ${BLAS_LIBRARIES}
    ${LAPACK_LIBRARIES}
)

# tuvx library
add_library(tuvx $<TARGET_OBJECTS:tuvx_object>)
add_library(musica::tuvx ALIAS tuvx)

set_target_properties(tuvx PROPERTIES 
  ARCHIVE_OUTPUT_DIRECTORY ${TUVX_LIB_DIR}
  VERSION ${PROJECT_VERSION}
  SOVERSION ${PROJECT_VERSION_MAJOR}
)

target_link_libraries(tuvx
  PUBLIC
    PkgConfig::netcdff
    yaml-cpp::yaml-cpp
    ${BLAS_LIBRARIES}
    ${LAPACK_LIBRARIES}
)

target_include_directories(tuvx
  PUBLIC  
    $<BUILD_INTERFACE:${TUVX_MOD_DIR}>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/../include>
    $<INSTALL_INTERFACE:${INSTALL_MOD_DIR}>
    $<INSTALL_INTERFACE:${INSTALL_INCLUDE_DIR}>
)

target_sources(tuvx_object
  PRIVATE
    constants.F90
    core.F90
    cross_section.F90
    cross_section_factory.F90
    cross_section_warehouse.F90
    diagnostic_util.F90
    dose_rates.F90
    grid.F90
    grid_factory.F90
    grid_warehouse.F90
    heating_rates.F90
    interpolate.F90
    la_sr_bands.F90
    linear_algebra.F90
    netcdf.F90
    output.F90
    photolysis_rates.F90
    profile.F90
    profile_factory.F90
    profile_warehouse.F90
    quantum_yield.F90
    quantum_yield_factory.F90
    spectral_weight.F90
    spectral_weight_factory.F90
    spherical_geometry.F90
    util.F90
)

add_subdirectory(cross_sections)
add_subdirectory(grids)
add_subdirectory(linear_algebras)
add_subdirectory(profiles)
add_subdirectory(quantum_yields)
add_subdirectory(radiative_transfer)
add_subdirectory(spectral_weights)
add_subdirectory(util)

################################################################################
