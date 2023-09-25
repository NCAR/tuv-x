include(CMakePackageConfigHelpers)

# install the binary and static library
install(
  TARGETS 
    tuv-x tuvx
  EXPORT 
    tuvx_Exports
  LIBRARY DESTINATION ${INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

# install the mod files
install(
  DIRECTORY 
    ${TUVX_MOD_DIR}/
  DESTINATION 
    ${INSTALL_MOD_DIR}
  FILES_MATCHING PATTERN "*.mod"
)

# install the data files
install(
  DIRECTORY 
    ${PROJECT_SOURCE_DIR}/data/
  DESTINATION 
    ${INSTALL_PREFIX}/data
)

# install the example config files
install(
  DIRECTORY 
    ${PROJECT_SOURCE_DIR}/examples/
  DESTINATION 
    ${INSTALL_PREFIX}/examples
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

# https://gitlab.kitware.com/cmake/community/-/wikis/FAQ#can-i-do-make-uninstall-with-cmake
if(NOT TARGET uninstall)
  configure_file(
    "${PROJECT_SOURCE_DIR}/cmake/cmake_uninstall.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake"
    IMMEDIATE @ONLY)

  add_custom_target(uninstall
    COMMAND ${CMAKE_COMMAND} -P ${CMAKE_CURRENT_BINARY_DIR}/cmake_uninstall.cmake)
endif()