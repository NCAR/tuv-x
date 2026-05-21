################################################################################
# TUV-x C++ library

add_library(tuvx INTERFACE)
add_library(musica::tuvx ALIAS tuvx)

target_compile_features(tuvx INTERFACE cxx_std_20)

# MSVC's <cmath> requires _USE_MATH_DEFINES to expose M_PI etc.
# Define it as a compile definition rather than in source to avoid the reserved-identifier warning.
if(MSVC)
  target_compile_definitions(tuvx INTERFACE _USE_MATH_DEFINES)
endif()

target_include_directories(
  tuvx
  INTERFACE
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

################################################################################
