# FindNetCDF.cmake
# Find the NetCDF C and Fortran libraries using nc-config/nf-config,
# with a pkg-config fallback.
#
# Imported targets created:
#   NetCDF::NetCDF_C       - NetCDF C library
#   NetCDF::NetCDF_Fortran - NetCDF Fortran library (links NetCDF::NetCDF_C)
#
# Result variables:
#   NetCDF_FOUND   - True if the requested components were found
#   NetCDF_VERSION - Version string of netcdf-c

include(FindPackageHandleStandardArgs)

# If a parent project already created the targets, accept them as-is
if(TARGET NetCDF::NetCDF_C AND TARGET NetCDF::NetCDF_Fortran)
  set(NetCDF_FOUND TRUE)
  return()
endif()

# ---- Locate via nc-config / nf-config (preferred; immune to PKG_CONFIG_PATH issues) ----

find_program(NC_CONFIG NAMES nc-config DOC "NetCDF C config script")
find_program(NF_CONFIG NAMES nf-config DOC "NetCDF Fortran config script")

if(NC_CONFIG)
  execute_process(COMMAND ${NC_CONFIG} --version
    OUTPUT_VARIABLE _nc_version OUTPUT_STRIP_TRAILING_WHITESPACE)
  string(REGEX REPLACE ".*([0-9]+\\.[0-9]+\\.[0-9]+).*" "\\1" NetCDF_VERSION "${_nc_version}")

  execute_process(COMMAND ${NC_CONFIG} --includedir
    OUTPUT_VARIABLE _nc_includedir OUTPUT_STRIP_TRAILING_WHITESPACE)
  execute_process(COMMAND ${NC_CONFIG} --libdir
    OUTPUT_VARIABLE _nc_libdir OUTPUT_STRIP_TRAILING_WHITESPACE)

  find_path(NetCDF_C_INCLUDE_DIR NAMES netcdf.h
    HINTS ${_nc_includedir} NO_DEFAULT_PATH)
  find_library(NetCDF_C_LIBRARY NAMES netcdf
    HINTS ${_nc_libdir} NO_DEFAULT_PATH)
endif()

if(NF_CONFIG)
  execute_process(COMMAND ${NF_CONFIG} --includedir
    OUTPUT_VARIABLE _nf_includedir OUTPUT_STRIP_TRAILING_WHITESPACE)
  execute_process(COMMAND ${NF_CONFIG} --libdir
    OUTPUT_VARIABLE _nf_libdir OUTPUT_STRIP_TRAILING_WHITESPACE)

  find_path(NetCDF_Fortran_INCLUDE_DIR NAMES netcdf.mod
    HINTS ${_nf_includedir} NO_DEFAULT_PATH)
  # nf-config --includedir is authoritative; use it directly if find_path missed the .mod file
  if(NOT NetCDF_Fortran_INCLUDE_DIR AND _nf_includedir)
    set(NetCDF_Fortran_INCLUDE_DIR "${_nf_includedir}"
      CACHE PATH "NetCDF Fortran include directory" FORCE)
  endif()

  find_library(NetCDF_Fortran_LIBRARY NAMES netcdff
    HINTS ${_nf_libdir} NO_DEFAULT_PATH)
endif()

# ---- pkg-config fallback ----

if(NOT NC_CONFIG OR NOT NetCDF_C_LIBRARY)
  find_package(PkgConfig QUIET)
  if(PkgConfig_FOUND)
    pkg_check_modules(_netcdf QUIET netcdf)
    if(_netcdf_FOUND)
      find_library(NetCDF_C_LIBRARY NAMES netcdf HINTS ${_netcdf_LIBRARY_DIRS})
      find_path(NetCDF_C_INCLUDE_DIR NAMES netcdf.h HINTS ${_netcdf_INCLUDE_DIRS})
    endif()
    pkg_check_modules(_netcdff QUIET netcdf-fortran)
    if(_netcdff_FOUND)
      find_library(NetCDF_Fortran_LIBRARY NAMES netcdff HINTS ${_netcdff_LIBRARY_DIRS})
      find_path(NetCDF_Fortran_INCLUDE_DIR NAMES netcdf.mod HINTS ${_netcdff_INCLUDE_DIRS})
      if(NOT NetCDF_Fortran_INCLUDE_DIR AND _netcdff_INCLUDE_DIRS)
        list(GET _netcdff_INCLUDE_DIRS 0 _first_incdir)
        set(NetCDF_Fortran_INCLUDE_DIR "${_first_incdir}"
          CACHE PATH "NetCDF Fortran include directory" FORCE)
      endif()
    endif()
  endif()
endif()

# ---- Determine required vars based on requested components ----

if(NOT NetCDF_FIND_COMPONENTS)
  set(NetCDF_FIND_COMPONENTS C Fortran)
endif()

set(_netcdf_required_vars NetCDF_C_LIBRARY NetCDF_C_INCLUDE_DIR)
foreach(_comp IN LISTS NetCDF_FIND_COMPONENTS)
  if(_comp STREQUAL "Fortran")
    list(APPEND _netcdf_required_vars NetCDF_Fortran_LIBRARY NetCDF_Fortran_INCLUDE_DIR)
  endif()
endforeach()

find_package_handle_standard_args(NetCDF
  REQUIRED_VARS ${_netcdf_required_vars}
  VERSION_VAR   NetCDF_VERSION
)

# ---- Create imported targets ----

if(NetCDF_FOUND)
  if(NOT TARGET NetCDF::NetCDF_C)
    add_library(NetCDF::NetCDF_C UNKNOWN IMPORTED)
    set_target_properties(NetCDF::NetCDF_C PROPERTIES
      IMPORTED_LOCATION             "${NetCDF_C_LIBRARY}"
      INTERFACE_INCLUDE_DIRECTORIES "${NetCDF_C_INCLUDE_DIR}"
    )
  endif()

  if(NetCDF_Fortran_LIBRARY AND NetCDF_Fortran_INCLUDE_DIR
      AND NOT TARGET NetCDF::NetCDF_Fortran)
    add_library(NetCDF::NetCDF_Fortran UNKNOWN IMPORTED)
    set_target_properties(NetCDF::NetCDF_Fortran PROPERTIES
      IMPORTED_LOCATION             "${NetCDF_Fortran_LIBRARY}"
      INTERFACE_INCLUDE_DIRECTORIES "${NetCDF_Fortran_INCLUDE_DIR}"
      INTERFACE_LINK_LIBRARIES      "NetCDF::NetCDF_C"
    )
  endif()
endif()

mark_as_advanced(
  NC_CONFIG NF_CONFIG
  NetCDF_C_INCLUDE_DIR NetCDF_C_LIBRARY
  NetCDF_Fortran_INCLUDE_DIR NetCDF_Fortran_LIBRARY
)
