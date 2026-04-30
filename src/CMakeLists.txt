################################################################################
# TUV-x C++ library

add_library(tuvx INTERFACE)
add_library(musica::tuvx ALIAS tuvx)

target_compile_features(tuvx INTERFACE cxx_std_20)

target_include_directories(
  tuvx
  INTERFACE
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)

################################################################################
