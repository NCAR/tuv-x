################################################################################
# Preamble

cmake_minimum_required(VERSION 3.17)

project(
  tuv-x
  VERSION 0.9.0
  LANGUAGES Fortran CXX C
)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH};${PROJECT_SOURCE_DIR}/cmake")

set(CMAKE_USER_MAKE_RULES_OVERRIDE "${PROJECT_SOURCE_DIR}/cmake/set_defaults.cmake")

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING
      "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel."
      FORCE)
endif(NOT CMAKE_BUILD_TYPE)

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
option(TUVX_ENABLE_MPI "Enable MPI parallel support" OFF)
cmake_dependent_option(TUVX_ENABLE_OPENMP "Enable OpenMP support" OFF "TUVX_ENABLE_MPI" OFF)
option(TUVX_ENABLE_LAPACK "Enable LAPACK" OFF)
option(TUVX_ENABLE_TESTS "Build tests" ON)
option(TUVX_ENABLE_BENCHMARK "Build benchmark examples" OFF)
option(TUVX_ENABLE_COVERAGE "Enable code coverage output" OFF)
option(TUVX_ENABLE_MEMCHECK "Enable memory checking in tests" OFF)
option(TUVX_ENABLE_NC_CONFIG "Use nc-config to determine NetCDF libraries" OFF)
option(TUVX_ENABLE_REGRESSION_TESTS "Enable regression tests" OFF)
option(TUVX_BUILD_DOCS "Build the documentation" OFF)
set(TUVX_MOD_DIR "${PROJECT_BINARY_DIR}/include" CACHE PATH "Directory to find Fortran module files during the build")
set(TUVX_INSTALL_MOD_DIR "${INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}" CACHE PATH "Directory to install Fortran module files")
set(TUVX_INSTALL_INCLUDE_DIR "${INSTALL_PREFIX}/${CMAKE_INSTALL_INCLUDEDIR}" CACHE PATH "Directory to install include files")

# Set up include and lib directories
set(TUVX_LIB_DIR "${PROJECT_BINARY_DIR}/lib")
include(GNUInstallDirs)
set(INSTALL_PREFIX "tuvx-${PROJECT_VERSION}" )

# Compiler specific flags
if(${CMAKE_Fortran_COMPILER_ID} MATCHES "NAG")
  add_definitions(-DMUSICA_IS_NAG_COMPILER)
elseif(${CMAKE_Fortran_COMPILER_ID} MATCHES "GNU")
  add_definitions(-DMUSICA_IS_GNU_COMPILER)
elseif(${CMAKE_Fortran_COMPILER_ID} MATCHES "Intel")
  add_definitions(-DMUSICA_IS_INTEL_COMPILER)
endif()

# MPI
if(TUVX_ENABLE_MPI)
  add_definitions(-DMUSICA_USE_MPI)
endif()

# copy data
if (TUVX_ENABLE_TESTS)
  add_custom_target(copy-data ALL COMMAND ${CMAKE_COMMAND}
    -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/data ${CMAKE_BINARY_DIR}/data)
  add_custom_target(copy-tool ALL COMMAND ${CMAKE_COMMAND}
    -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/tool ${CMAKE_BINARY_DIR}/tool)
  add_custom_target(copy-examples ALL COMMAND ${CMAKE_COMMAND}
    -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/examples ${CMAKE_BINARY_DIR}/examples)
endif()

################################################################################
# Dependencies

include(cmake/dependencies.cmake)

################################################################################
# TUV-x targets

add_subdirectory(src)

# add tuv-x standalone executable
configure_file(src/version.F90.in version.F90 @ONLY)

add_executable(tuv-x src/tuvx.F90 version.F90)

set_target_properties(tuv-x
  PROPERTIES
  LINKER_LANGUAGE Fortran
)

if(NOT BUILD_SHARED_LIBS)
  set_target_properties(tuv-x
    PROPERTIES
    POSITION_INDEPENDENT_CODE ON
  )
endif()

target_link_libraries(tuv-x 
  PUBLIC 
    musica::tuvx
    yaml-cpp::yaml-cpp
    ${BLAS_LIBRARIES}
)

target_include_directories(tuv-x
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}/include>
    $<INSTALL_INTERFACE:${TUVX_INSTALL_INCLUDE_DIR}>)

if(TUVX_ENABLE_OPENMP)
  target_link_libraries(tuv-x PUBLIC OpenMP::OpenMP_Fortran)
endif()

if(TUVX_ENABLE_LAPACK)
  target_link_libraries(tuv-x PUBLIC ${LAPACK_LIBRARIES})
endif()
#   set(${LAPACK_LIBRARIES} "")
# endif()

################################################################################
# TUV-x docs

if(TUVX_BUILD_DOCS)
  add_subdirectory(docs)
endif()

################################################################################
# TUV-x tests

if(PROJECT_IS_TOP_LEVEL AND TUVX_ENABLE_TESTS)
  # Test code coverage
  if(TUVX_ENABLE_COVERAGE)
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
# benchmarking

if(TUVX_ENABLE_BENCHMARK AND TUVX_ENABLE_LAPACK)
  add_subdirectory(benchmark)
endif()

################################################################################
