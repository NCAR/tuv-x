find_package(PkgConfig REQUIRED)
include(FetchContent)

################################################################################
# LAPACK

find_package(BLAS)
find_package(LAPACK)

################################################################################
# Memory check

if(ENABLE_MEMCHECK)
  find_file(MEMCHECK_SUPPRESS_FILE
    DOC "Suppression file for memory checking"
    NAMES openmpi-valgrind.supp
    PATHS
      /usr/share/openmpi
      /usr/lib64/openmpi/share
      /usr/lib64/openmpi/share/openmpi
      /usr/share)
  if(MEMCHECK_SUPPRESS_FILE)
    set(MEMCHECK_SUPPRESS "--suppressions=${PROJECT_SOURCE_DIR}/test/valgrind.supp --suppressions=${MEMCHECK_SUPPRESS_FILE}")
  else()
    set(MEMCHECK_SUPPRESS "--suppressions=${PROJECT_SOURCE_DIR}/test/valgrind.supp")
  endif()
endif()

################################################################################
# OpenMP

if(ENABLE_OPENMP)
  find_package(OpenMP)
  if(OpenMP_Fortran_FOUND)
    message(STATUS "Compiling with OpenMP support")
    add_definitions(-DMUSICA_USE_OPENMP)
  else()
    message(FATAL_ERROR "OpenMP package not found")
  endif()
endif()

################################################################################
# NetCDF library

find_package(PkgConfig REQUIRED)
pkg_check_modules(netcdff IMPORTED_TARGET REQUIRED netcdf-fortran)

################################################################################
# yaml-cpp

FetchContent_Declare(
  yaml-cpp
  GIT_REPOSITORY https://github.com/jbeder/yaml-cpp/
  GIT_TAG 0.8.0
)
FetchContent_MakeAvailable(yaml-cpp)

################################################################################
# Docs

if(BUILD_DOCS)
  find_package(Doxygen REQUIRED)
  find_package(Sphinx REQUIRED)
endif()
