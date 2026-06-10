include(FetchContent)

# ##############################################################################
# LAPACK

if(TUVX_ENABLE_LAPACK AND NOT TUVX_DOCS_ONLY)
  find_package(LAPACK)
  find_package(LAPACKE)
  find_package(BLAS)
endif()

# ##############################################################################
# Memory check

if(TUVX_ENABLE_MEMCHECK AND NOT TUVX_DOCS_ONLY)
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

if(TUVX_ENABLE_OPENMP AND NOT TUVX_DOCS_ONLY)
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
#
# When TUV-x is built inside a host project (e.g. MUSICA, or a model embedding
# MUSICA such as CATChem) NetCDF has usually already been located and exposed as
# an imported target. Different finders use different target names, so probe the
# common variants and reuse a host-provided target before searching ourselves.
# Only when nothing is found do we fall back to pkg-config. The resolved target
# names are stored in TUVX_NETCDF_C_TARGET / TUVX_NETCDF_FORTRAN_TARGET and linked
# directly (we avoid wrapping them in an imported target, which would leak an
# unresolvable target into the exported link interface).

if(NOT TUVX_DOCS_ONLY)
  # Return the first of ARGN that already exists as a target (target names are
  # case-sensitive, so we list each convention explicitly).
  function(tuvx_first_existing_target out_var)
    foreach(candidate ${ARGN})
      if(TARGET ${candidate})
        set(${out_var} ${candidate} PARENT_SCOPE)
        return()
      endif()
    endforeach()
    set(${out_var} "" PARENT_SCOPE)
  endfunction()

  tuvx_first_existing_target(TUVX_NETCDF_C_TARGET
    netCDF::netcdf
    NetCDF::NetCDF_C
    netcdf
    NETCDF::NETCDF)
  tuvx_first_existing_target(TUVX_NETCDF_FORTRAN_TARGET
    netCDF::netcdff
    NetCDF::NetCDF_Fortran
    netcdff)

  if(TUVX_NETCDF_C_TARGET AND TUVX_NETCDF_FORTRAN_TARGET)
    message(STATUS "TUV-x: reusing NetCDF targets from parent project: "
                   "${TUVX_NETCDF_C_TARGET}, ${TUVX_NETCDF_FORTRAN_TARGET}")
  else()
    message(STATUS "TUV-x: no parent NetCDF target found; locating via pkg-config")
    find_package(PkgConfig REQUIRED)
    pkg_check_modules(netcdff IMPORTED_TARGET REQUIRED netcdf-fortran)
    pkg_check_modules(netcdfc IMPORTED_TARGET REQUIRED netcdf)
    set(TUVX_NETCDF_C_TARGET PkgConfig::netcdfc)
    set(TUVX_NETCDF_FORTRAN_TARGET PkgConfig::netcdff)
  endif()
endif()

# ##############################################################################
# yaml-cpp
if(NOT TUVX_DOCS_ONLY)
  FetchContent_Declare(
    yaml-cpp
    GIT_REPOSITORY https://github.com/jbeder/yaml-cpp/
    GIT_TAG 65c1c270dbe7eec37b2df2531d7497c4eea79aee
    GIT_PROGRESS NOT
    FIND_PACKAGE_ARGS NAMES yaml-cpp
    ${FETCHCONTENT_QUIET})

  set(YAML_CPP_BUILD_TOOLS OFF CACHE BOOL "" FORCE)

  FetchContent_MakeAvailable(yaml-cpp)

  # Ensure yaml-cpp::yaml-cpp target exists (handle both system and FetchContent scenarios)
  if(NOT TARGET yaml-cpp::yaml-cpp)
    if(TARGET yaml-cpp)
      # Create alias for system-installed yaml-cpp that provides 'yaml-cpp' target
      add_library(yaml-cpp::yaml-cpp ALIAS yaml-cpp)
    endif()
  endif()
endif()

# ##############################################################################
# Docs

if(TUVX_BUILD_DOCS)
  find_package(Doxygen REQUIRED)
  find_package(Sphinx REQUIRED)
endif()

# ##############################################################################
# google test and benchmark

if(TUVX_ENABLE_BENCHMARK AND NOT TUVX_DOCS_ONLY)
  FetchContent_Declare(
    googlebenchmark
    GIT_REPOSITORY https://github.com/google/benchmark.git
    GIT_TAG v1.8.3
    FIND_PACKAGE_ARGS NAMES benchmark)

  set(BENCHMARK_DOWNLOAD_DEPENDENCIES ON)
  set(BENCHMARK_ENABLE_GTEST_TESTS OFF)
  set(BENCHMARK_ENABLE_ASSEMBLY_TESTS OFF)
  set(BENCHMARK_ENABLE_TESTING OFF)
  FetchContent_MakeAvailable(googlebenchmark)
  
  # Ensure benchmark::benchmark target exists (handle both system and FetchContent scenarios)
  if(NOT TARGET benchmark::benchmark)
    if(TARGET benchmark)
      # Create alias for system-installed benchmark that provides 'benchmark' target
      add_library(benchmark::benchmark ALIAS benchmark)
    endif()
  endif()
endif()

if(TUVX_ENABLE_TESTS AND NOT TUVX_DOCS_ONLY)
  FetchContent_Declare(
    googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG be03d00f5f0cc3a997d1a368bee8a1fe93651f48
    FIND_PACKAGE_ARGS NAMES GTest)

  set(INSTALL_GTEST
      OFF
      CACHE BOOL "" FORCE)
  set(BUILD_GMOCK
      OFF
      CACHE BOOL "" FORCE)

  FetchContent_MakeAvailable(googletest)

  # Ensure GTest::gtest_main target exists (handle both system and FetchContent scenarios)
  if(NOT TARGET GTest::gtest_main)
    if(TARGET gtest_main)
      # Create alias for system-installed gtest that provides 'gtest_main' target
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

# ##############################################################################
