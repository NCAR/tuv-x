include(FetchContent)

################################################################################
# LAPACK

if(TUVX_ENABLE_LAPACK)
  find_package(LAPACK)
  find_package(LAPACKE)
  find_package(BLAS)
endif()

################################################################################
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
        "--suppressions=${PROJECT_SOURCE_DIR}/test/valgrind.supp --suppressions=${MEMCHECK_SUPPRESS_FILE}")
  else()
    set(MEMCHECK_SUPPRESS
        "--suppressions=${PROJECT_SOURCE_DIR}/test/valgrind.supp")
  endif()
endif()

################################################################################
# NetCDF-C (optional — needed for data readers in Phase 5+)

if(TUVX_ENABLE_NETCDF)
  find_package(PkgConfig REQUIRED)
  pkg_check_modules(netcdfc IMPORTED_TARGET REQUIRED netcdf)
endif()

################################################################################
# Docs

if(TUVX_BUILD_DOCS)
  find_package(Doxygen REQUIRED)
  find_package(Sphinx REQUIRED)
endif()

################################################################################
# Google Benchmark

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

  if(NOT TARGET benchmark::benchmark)
    if(TARGET benchmark)
      add_library(benchmark::benchmark ALIAS benchmark)
    endif()
  endif()
endif()

################################################################################
# Google Test

if(TUVX_ENABLE_TESTS)
  FetchContent_Declare(
    googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG be03d00f5f0cc3a997d1a368bee8a1fe93651f48)

  set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)
  set(BUILD_GMOCK OFF CACHE BOOL "" FORCE)
  set(BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE)

  FetchContent_MakeAvailable(googletest)

  if(NOT TARGET GTest::gtest_main)
    if(TARGET gtest_main)
      add_library(GTest::gtest_main ALIAS gtest_main)
    endif()
  endif()

  # don't run clang-tidy on google test (only when built from source)
  if(TARGET gtest)
    set_target_properties(gtest PROPERTIES CXX_CLANG_TIDY "")
  endif()
  if(TARGET gtest_main)
    set_target_properties(gtest_main PROPERTIES CXX_CLANG_TIDY "")
  endif()
endif()

################################################################################
