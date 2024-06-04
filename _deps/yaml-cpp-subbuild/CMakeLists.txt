# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION 3.26.3)

# We name the project and the target for the ExternalProject_Add() call
# to something that will highlight to the user what we are working on if
# something goes wrong and an error message is produced.

project(yaml-cpp-populate NONE)


# Pass through things we've already detected in the main project to avoid
# paying the cost of redetecting them again in ExternalProject_Add()
set(GIT_EXECUTABLE [==[/glade/u/apps/casper/23.10/opt/view/bin/git]==])
set(GIT_VERSION_STRING [==[2.41.0]==])
set_property(GLOBAL PROPERTY _CMAKE_FindGit_GIT_EXECUTABLE_VERSION
  [==[/glade/u/apps/casper/23.10/opt/view/bin/git;2.41.0]==]
)


include(ExternalProject)
ExternalProject_Add(yaml-cpp-populate
                     "UPDATE_DISCONNECTED" "False" "GIT_REPOSITORY" "https://github.com/jbeder/yaml-cpp/" "GIT_TAG" "0.8.0" "GIT_PROGRESS" "NOT" "ON"
                    SOURCE_DIR          "/glade/derecho/scratch/adityad/tuv-x/_deps/yaml-cpp-src"
                    BINARY_DIR          "/glade/derecho/scratch/adityad/tuv-x/_deps/yaml-cpp-build"
                    CONFIGURE_COMMAND   ""
                    BUILD_COMMAND       ""
                    INSTALL_COMMAND     ""
                    TEST_COMMAND        ""
                    USES_TERMINAL_DOWNLOAD  YES
                    USES_TERMINAL_UPDATE    YES
                    USES_TERMINAL_PATCH     YES
)


