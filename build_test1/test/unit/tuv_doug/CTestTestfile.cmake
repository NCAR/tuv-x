# CMake generated Testfile for 
# Source directory: /home/runner/work/tuv-x/tuv-x/test/unit/tuv_doug
# Build directory: /home/runner/work/tuv-x/tuv-x/build_test1/test/unit/tuv_doug
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(data_sets "/home/runner/work/tuv-x/tuv-x/build_test1/test_data_sets")
set_tests_properties(data_sets PROPERTIES  WORKING_DIRECTORY "/home/runner/work/tuv-x/tuv-x/build_test1" _BACKTRACE_TRIPLES "/home/runner/work/tuv-x/tuv-x/cmake/test_util.cmake;83;add_test;/home/runner/work/tuv-x/tuv-x/test/unit/tuv_doug/CMakeLists.txt;43;add_tuvx_test;/home/runner/work/tuv-x/tuv-x/test/unit/tuv_doug/CMakeLists.txt;0;")
add_test(la_srb_lut "/home/runner/work/tuv-x/tuv-x/build_test1/test_la_srb_lut")
set_tests_properties(la_srb_lut PROPERTIES  WORKING_DIRECTORY "/home/runner/work/tuv-x/tuv-x/build_test1" _BACKTRACE_TRIPLES "/home/runner/work/tuv-x/tuv-x/cmake/test_util.cmake;83;add_test;/home/runner/work/tuv-x/tuv-x/test/unit/tuv_doug/CMakeLists.txt;55;add_tuvx_test;/home/runner/work/tuv-x/tuv-x/test/unit/tuv_doug/CMakeLists.txt;0;")
subdirs("JCALC")
