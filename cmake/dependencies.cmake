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
# musica-core library

if(${CMAKE_VERSION} VERSION_LESS "3.24") 
    find_package(musicacore REQUIRED)
else()
  include(FetchContent)

  set(ENABLE_UTIL_ONLY ON)

  FetchContent_Declare(musicacore
    GIT_REPOSITORY https://github.com/NCAR/musica-core.git
    GIT_TAG develop-fix-yaml-config # v0.4.2
    FIND_PACKAGE_ARGS NAMES musicacore
  )

  FetchContent_MakeAvailable(musicacore)
endif()

################################################################################
# Docs

if(BUILD_DOCS)
  find_package(Doxygen REQUIRED)
  find_package(Sphinx REQUIRED)
endif()
