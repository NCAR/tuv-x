# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION 3.5)

file(MAKE_DIRECTORY
  "/glade/derecho/scratch/adityad/tuv-x/_deps/yaml-cpp-src"
  "/glade/derecho/scratch/adityad/tuv-x/_deps/yaml-cpp-build"
  "/glade/derecho/scratch/adityad/tuv-x/_deps/yaml-cpp-subbuild/yaml-cpp-populate-prefix"
  "/glade/derecho/scratch/adityad/tuv-x/_deps/yaml-cpp-subbuild/yaml-cpp-populate-prefix/tmp"
  "/glade/derecho/scratch/adityad/tuv-x/_deps/yaml-cpp-subbuild/yaml-cpp-populate-prefix/src/yaml-cpp-populate-stamp"
  "/glade/derecho/scratch/adityad/tuv-x/_deps/yaml-cpp-subbuild/yaml-cpp-populate-prefix/src"
  "/glade/derecho/scratch/adityad/tuv-x/_deps/yaml-cpp-subbuild/yaml-cpp-populate-prefix/src/yaml-cpp-populate-stamp"
)

set(configSubDirs )
foreach(subDir IN LISTS configSubDirs)
    file(MAKE_DIRECTORY "/glade/derecho/scratch/adityad/tuv-x/_deps/yaml-cpp-subbuild/yaml-cpp-populate-prefix/src/yaml-cpp-populate-stamp/${subDir}")
endforeach()
if(cfgdir)
  file(MAKE_DIRECTORY "/glade/derecho/scratch/adityad/tuv-x/_deps/yaml-cpp-subbuild/yaml-cpp-populate-prefix/src/yaml-cpp-populate-stamp${cfgdir}") # cfgdir has leading slash
endif()
