include(CMakePackageConfigHelpers)

# install the header-only library
install(
  TARGETS
    tuvx
  EXPORT
    tuvx_Exports
)

# install headers
install(
  DIRECTORY
    ${PROJECT_SOURCE_DIR}/include/
  DESTINATION
    ${CMAKE_INSTALL_INCLUDEDIR}
)

# install the cmake config files
set(cmake_config_install_location ${INSTALL_PREFIX}/cmake/tuvx-${PROJECT_VERSION})

install(
  EXPORT
    tuvx_Exports
  DESTINATION
    ${cmake_config_install_location}
  NAMESPACE musica::
)

configure_package_config_file(
  "${PROJECT_SOURCE_DIR}/cmake/tuvxConfig.cmake.in"
  "${PROJECT_BINARY_DIR}/tuvxConfig.cmake"
  INSTALL_DESTINATION
    ${cmake_config_install_location}
)

write_basic_package_version_file(
  "${PROJECT_BINARY_DIR}/tuvxConfigVersion.cmake"
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY SameMajorVersion
)

install(
  FILES
    ${PROJECT_BINARY_DIR}/tuvxConfig.cmake
    ${PROJECT_BINARY_DIR}/tuvxConfigVersion.cmake
  DESTINATION
    ${INSTALL_PREFIX}/cmake/tuvx-${PROJECT_VERSION}
)

######################################################################
# uninstall target

if(NOT TARGET uninstall)
  configure_file(
    "${PROJECT_SOURCE_DIR}/cmake/cmake_uninstall.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
    IMMEDIATE @ONLY)

  add_custom_target(uninstall
    COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)
endif()
