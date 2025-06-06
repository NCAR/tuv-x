find_package(PkgConfig REQUIRED)
include(FetchContent)

# ##############################################################################
# LAPACK

if(TUVX_ENABLE_LAPACK)
  find_package(LAPACK)
  find_package(LAPACKE)
  find_package(BLAS)
endif()

# ##############################################################################
# Memory check

if(TUVX_ENABLE_MEMCHECK)
  find_file(
    MEMCHECK_SUPPRESS_FILE
    DOC "Suppression file for memory checking"
    NAMES openmpi-valgrind.supp
    PATHS /usr/share/openmpi /usr/lib64/openmpi/share
          /usr/lib64/openmpi/share/openmpi /usr/share)
  if(MEMCHECK_SUPPRESS_FILE)
    set(MEMCHECK_SUPPRESS
        "--suppressions=${PROJECT_SOURCE_DIR}/test/valgrind.supp --suppressions=${MEMCHECK_SUPPRESS_FILE}"
    )
  else()
    set(MEMCHECK_SUPPRESS
        "--suppressions=${PROJECT_SOURCE_DIR}/test/valgrind.supp")
  endif()
endif()

# ##############################################################################
# OpenMP

if(TUVX_ENABLE_OPENMP)
  find_package(OpenMP)
  if(OpenMP_Fortran_FOUND)
    message(STATUS "Compiling with OpenMP support")
    add_definitions(-DMUSICA_USE_OPENMP)
  else()
    message(FATAL_ERROR "OpenMP package not found")
  endif()
endif()

# ##############################################################################
# NetCDF library

find_package(PkgConfig REQUIRED)

pkg_check_modules(netcdff IMPORTED_TARGET REQUIRED netcdf-fortran)
pkg_check_modules(netcdfc IMPORTED_TARGET REQUIRED netcdf)

# ##############################################################################
# yaml-cpp

FetchContent_Declare(
  yaml-cpp
  GIT_REPOSITORY https://github.com/jbeder/yaml-cpp/
  GIT_TAG 28f93bdec6387d42332220afa9558060c8016795
  GIT_PROGRESS NOT
  ${FETCHCONTENT_QUIET})

set(YAML_CPP_BUILD_TOOLS OFF CACHE BOOL "" FORCE)

FetchContent_MakeAvailable(yaml-cpp)

# ##############################################################################
# Docs

if(TUVX_BUILD_DOCS)
  find_package(Doxygen REQUIRED)
  find_package(Sphinx REQUIRED)
endif()

# ##############################################################################
# google test and benchmark

if(TUVX_ENABLE_BENCHMARK)
  FetchContent_Declare(
    googlebenchmark
    GIT_REPOSITORY https://github.com/google/benchmark.git
    GIT_TAG v1.8.3)

  set(BENCHMARK_DOWNLOAD_DEPENDENCIES ON)
  set(BENCHMARK_ENABLE_GTEST_TESTS OFF)
  set(BENCHMARK_ENABLE_ASSEMBLY_TESTS OFF)
  set(BENCHMARK_ENABLE_TESTING OFF)
  FetchContent_MakeAvailable(googlebenchmark)
endif()

if(TUVX_ENABLE_TESTS)
  FetchContent_Declare(
    googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG be03d00f5f0cc3a997d1a368bee8a1fe93651f48)

  set(INSTALL_GTEST
      OFF
      CACHE BOOL "" FORCE)
  set(BUILD_GMOCK
      OFF
      CACHE BOOL "" FORCE)

  FetchContent_MakeAvailable(googletest)

  # don't run clang-tidy on google test
  set_target_properties(gtest PROPERTIES CXX_CLANG_TIDY "")
  set_target_properties(gtest_main PROPERTIES CXX_CLANG_TIDY "")
endif()

# ##############################################################################
