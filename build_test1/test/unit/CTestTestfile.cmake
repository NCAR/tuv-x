# CMake generated Testfile for 
# Source directory: /home/runner/work/tuv-x/tuv-x/test/unit
# Build directory: /home/runner/work/tuv-x/tuv-x/build_test1/test/unit
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(cxx_grid "/home/runner/work/tuv-x/tuv-x/build_test1/test_cxx_grid")
set_tests_properties(cxx_grid PROPERTIES  WORKING_DIRECTORY "/home/runner/work/tuv-x/tuv-x/build_test1" _BACKTRACE_TRIPLES "/home/runner/work/tuv-x/tuv-x/cmake/test_util.cmake;83;add_test;/home/runner/work/tuv-x/tuv-x/cmake/test_util.cmake;71;add_tuvx_test;/home/runner/work/tuv-x/tuv-x/test/unit/CMakeLists.txt;23;create_standard_cxx_test;/home/runner/work/tuv-x/tuv-x/test/unit/CMakeLists.txt;0;")
add_test(cxx_profile "/home/runner/work/tuv-x/tuv-x/build_test1/test_cxx_profile")
set_tests_properties(cxx_profile PROPERTIES  WORKING_DIRECTORY "/home/runner/work/tuv-x/tuv-x/build_test1" _BACKTRACE_TRIPLES "/home/runner/work/tuv-x/tuv-x/cmake/test_util.cmake;83;add_test;/home/runner/work/tuv-x/tuv-x/cmake/test_util.cmake;71;add_tuvx_test;/home/runner/work/tuv-x/tuv-x/test/unit/CMakeLists.txt;24;create_standard_cxx_test;/home/runner/work/tuv-x/tuv-x/test/unit/CMakeLists.txt;0;")
add_test(grid_warehouse "/home/runner/work/tuv-x/tuv-x/build_test1/test_grid_warehouse")
set_tests_properties(grid_warehouse PROPERTIES  WORKING_DIRECTORY "/home/runner/work/tuv-x/tuv-x/build_test1" _BACKTRACE_TRIPLES "/home/runner/work/tuv-x/tuv-x/cmake/test_util.cmake;83;add_test;/home/runner/work/tuv-x/tuv-x/cmake/test_util.cmake;39;add_tuvx_test;/home/runner/work/tuv-x/tuv-x/test/unit/CMakeLists.txt;25;create_standard_test;/home/runner/work/tuv-x/tuv-x/test/unit/CMakeLists.txt;0;")
add_test(heating_rates "/home/runner/work/tuv-x/tuv-x/build_test1/test_heating_rates")
set_tests_properties(heating_rates PROPERTIES  WORKING_DIRECTORY "/home/runner/work/tuv-x/tuv-x/build_test1" _BACKTRACE_TRIPLES "/home/runner/work/tuv-x/tuv-x/cmake/test_util.cmake;83;add_test;/home/runner/work/tuv-x/tuv-x/cmake/test_util.cmake;39;add_tuvx_test;/home/runner/work/tuv-x/tuv-x/test/unit/CMakeLists.txt;26;create_standard_test;/home/runner/work/tuv-x/tuv-x/test/unit/CMakeLists.txt;0;")
add_test(la_sr_bands "/home/runner/work/tuv-x/tuv-x/build_test1/test_la_sr_bands")
set_tests_properties(la_sr_bands PROPERTIES  WORKING_DIRECTORY "/home/runner/work/tuv-x/tuv-x/build_test1" _BACKTRACE_TRIPLES "/home/runner/work/tuv-x/tuv-x/cmake/test_util.cmake;83;add_test;/home/runner/work/tuv-x/tuv-x/cmake/test_util.cmake;39;add_tuvx_test;/home/runner/work/tuv-x/tuv-x/test/unit/CMakeLists.txt;27;create_standard_test;/home/runner/work/tuv-x/tuv-x/test/unit/CMakeLists.txt;0;")
add_test(spherical_geometry "/home/runner/work/tuv-x/tuv-x/build_test1/test_spherical_geometry")
set_tests_properties(spherical_geometry PROPERTIES  WORKING_DIRECTORY "/home/runner/work/tuv-x/tuv-x/build_test1" _BACKTRACE_TRIPLES "/home/runner/work/tuv-x/tuv-x/cmake/test_util.cmake;83;add_test;/home/runner/work/tuv-x/tuv-x/cmake/test_util.cmake;39;add_tuvx_test;/home/runner/work/tuv-x/tuv-x/test/unit/CMakeLists.txt;28;create_standard_test;/home/runner/work/tuv-x/tuv-x/test/unit/CMakeLists.txt;0;")
subdirs("cross_section")
subdirs("grid")
subdirs("profile")
subdirs("quantum_yield")
subdirs("radiative_transfer")
subdirs("radiator")
subdirs("spectral_weight")
subdirs("tuv_doug")
subdirs("linear_algebra")
subdirs("util")
