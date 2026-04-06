################################################################################
# Preamble

cmake_minimum_required(VERSION 3.21)

# project and version must be on the same line so that the docs can extract it
project(tuv-x VERSION 0.15.0 LANGUAGES CXX C)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})
set(CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH};${PROJECT_SOURCE_DIR}/cmake")

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING
      "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel."
      FORCE)
endif()

message(STATUS "CMake build configuration for tuv-x(${CMAKE_BUILD_TYPE}) ${PROJECT_VERSION}")

################################################################################
# Project wide setup options

# Placeholders — no enable_language() or target logic yet. Wired in when DeviceArrayPolicy lands.
option(TUVX_ENABLE_CUDA      "Enable CUDA support"           OFF)
option(TUVX_ENABLE_HIP       "Enable HIP support"            OFF)
option(TUVX_ENABLE_TESTS     "Build unit tests"              ON)
option(TUVX_ENABLE_BENCHMARK "Build benchmarks"              OFF)
option(TUVX_ENABLE_COVERAGE  "Enable code coverage output"   OFF)
option(TUVX_ENABLE_MEMCHECK  "Enable Valgrind memcheck"      OFF)
option(TUVX_ENABLE_LAPACK    "Enable LAPACK for benchmarks"  OFF)
option(TUVX_ENABLE_NETCDF    "Enable NetCDF-C support"       OFF)
option(TUVX_BUILD_DOCS       "Build the documentation"       OFF)

# Set up include and lib directories
set(TUVX_LIB_DIR "${PROJECT_BINARY_DIR}/lib")
include(GNUInstallDirs)
set(INSTALL_PREFIX "tuvx-${PROJECT_VERSION}")

################################################################################
# Dependencies

include(cmake/dependencies.cmake)

################################################################################
# TUV-x library target

add_subdirectory(src)

################################################################################
# TUV-x docs

if(TUVX_BUILD_DOCS)
  add_subdirectory(docs)
endif()

################################################################################
# TUV-x tests

if(PROJECT_IS_TOP_LEVEL AND TUVX_ENABLE_TESTS)
  if(TUVX_ENABLE_COVERAGE)
    include(CodeCoverage)
    setup_target_for_coverage_lcov(
        NAME coverage
        EXECUTABLE "ctest"
        EXCLUDE "${PROJECT_SOURCE_DIR}/test/*"
        BASE_DIRECTORY "${PROJECT_SOURCE_DIR}/src")
  endif()

  enable_testing()
  add_subdirectory(test)
endif()

################################################################################
# Packaging

if(PROJECT_IS_TOP_LEVEL)
  add_subdirectory(packaging)
endif()

################################################################################
# Benchmarking

if(TUVX_ENABLE_BENCHMARK AND TUVX_ENABLE_LAPACK)
  add_subdirectory(benchmark)
endif()

################################################################################
