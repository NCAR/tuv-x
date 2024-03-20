################################################################################
# Preamble

cmake_minimum_required(VERSION 3.17)
set(CMAKE_USER_MAKE_RULES_OVERRIDE "cmake/set_defaults.cmake")

project(
  tuv-x
  VERSION 0.8.0
  LANGUAGES Fortran CXX C
)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH};${PROJECT_SOURCE_DIR}/cmake")

message ( STATUS "CMake build configuration for tuv-x(${CMAKE_BUILD_TYPE}) ${PROJECT_VERSION}" )

if(${CMAKE_VERSION} VERSION_LESS "3.21") 
  if(CMAKE_CURRENT_SOURCE_DIR STREQUAL CMAKE_SOURCE_DIR)
    set(PROJECT_IS_TOP_LEVEL TRUE)
  else()
    set(PROJECT_IS_TOP_LEVEL FALSE)
  endif()
endif()

################################################################################
# Projet wide setup options

include(CMakeDependentOption)
option(ENABLE_MPI "Enable MPI parallel support" OFF)
cmake_dependent_option(ENABLE_OPENMP "Enable OpenMP support" OFF "ENABLE_MPI" OFF)
option(ENABLE_TESTS "Build tests" ON)
option(ENABLE_COVERAGE "Enable code coverage output" OFF)
option(ENABLE_MEMCHECK "Enable memory checking in tests" OFF)
option(ENABLE_NC_CONFIG "Use nc-config to determine NetCDF libraries" OFF)
option(BUILD_DOCS "Build the documentation" OFF)

# Set up include and lib directories
set(TUVX_MOD_DIR "${PROJECT_BINARY_DIR}/include")
set(TUVX_LIB_DIR "${PROJECT_BINARY_DIR}/lib")
include(GNUInstallDirs)
set(INSTALL_PREFIX "tuvx-${PROJECT_VERSION}" )
set(INSTALL_MOD_DIR "${INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}")

# Compiler specific flags
if(${CMAKE_Fortran_COMPILER_ID} MATCHES "NAG")
  add_definitions(-DMUSICA_IS_NAG_COMPILER)
endif()

# MPI
if(ENABLE_MPI)
  add_definitions(-DMUSICA_USE_MPI)
endif()

# copy data
add_custom_target(copy-data ALL COMMAND ${CMAKE_COMMAND}
  -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/data ${CMAKE_BINARY_DIR}/data)
add_custom_target(copy-tool ALL COMMAND ${CMAKE_COMMAND}
  -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/tool ${CMAKE_BINARY_DIR}/tool)
add_custom_target(copy-examples ALL COMMAND ${CMAKE_COMMAND}
  -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/examples ${CMAKE_BINARY_DIR}/examples)

################################################################################
# Dependencies

include(cmake/dependencies.cmake)

################################################################################
# TUV-x targets

add_subdirectory(src)

# add tuv-x standalone executable
configure_file(src/version.F90.in version.F90 @ONLY)

add_executable(tuv-x src/tuvx.F90 version.F90)

target_link_libraries(tuv-x 
  PUBLIC 
    musica::tuvx
    yaml-cpp::yaml-cpp
    ${BLAS_LIBRARIES}
    ${LAPACK_LIBRARIES}
)

target_include_directories(tuv-x
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/include>
    $<INSTALL_INTERFACE:${INSTALL_INCLUDE_DIR}>)

if(ENABLE_OPENMP)
  target_link_libraries(tuv-x PUBLIC OpenMP::OpenMP_Fortran)
endif()

################################################################################
# TUV-x docs

if(BUILD_DOCS)
  add_subdirectory(docs)
endif()

################################################################################
# TUV-x tests

if(PROJECT_IS_TOP_LEVEL AND ENABLE_TESTS)
  # Test code coverage
  if(ENABLE_COVERAGE)
    include(CodeCoverage)
    if(${CMAKE_Fortran_COMPILER_ID} MATCHES "GNU")
      setup_target_for_coverage_lcov(
          NAME coverage
          EXECUTABLE "ctest"
          EXCLUDE "${PROJECT_SOURCE_DIR}/test/*"
          BASE_DIRECTORY "${PROJECT_SOURCE_DIR}/src")
    endif()
  endif()

  enable_testing()
  add_subdirectory(test)
endif()

################################################################################
# Packaging

# only include packaging if we are the top level project being built
if(PROJECT_IS_TOP_LEVEL)
  add_subdirectory(packaging)
endif()

################################################################################