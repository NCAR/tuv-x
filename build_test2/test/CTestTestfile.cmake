# CMake generated Testfile for 
# Source directory: /home/runner/work/tuv-x/tuv-x/test
# Build directory: /home/runner/work/tuv-x/tuv-x/build_test2/test
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(TUV_5_4 "/home/runner/work/tuv-x/tuv-x/build_test2/tuv-x" "../examples/tuv_5_4.json")
set_tests_properties(TUV_5_4 PROPERTIES  WORKING_DIRECTORY "/home/runner/work/tuv-x/tuv-x/build_test2/example_tuv_5_4" _BACKTRACE_TRIPLES "/home/runner/work/tuv-x/tuv-x/test/CMakeLists.txt;47;add_test;/home/runner/work/tuv-x/tuv-x/test/CMakeLists.txt;0;")
add_test(TS1_TSMLT "/home/runner/work/tuv-x/tuv-x/build_test2/tuv-x" "../examples/ts1_tsmlt.json")
set_tests_properties(TS1_TSMLT PROPERTIES  WORKING_DIRECTORY "/home/runner/work/tuv-x/tuv-x/build_test2/example_ts1_tsmlt" _BACKTRACE_TRIPLES "/home/runner/work/tuv-x/tuv-x/test/CMakeLists.txt;53;add_test;/home/runner/work/tuv-x/tuv-x/test/CMakeLists.txt;0;")
add_test(TUV_5_4_YAML "/home/runner/work/tuv-x/tuv-x/build_test2/tuv-x" "../examples/tuv_5_4.yml")
set_tests_properties(TUV_5_4_YAML PROPERTIES  WORKING_DIRECTORY "/home/runner/work/tuv-x/tuv-x/build_test2/example_tuv_5_4_yaml" _BACKTRACE_TRIPLES "/home/runner/work/tuv-x/tuv-x/test/CMakeLists.txt;60;add_test;/home/runner/work/tuv-x/tuv-x/test/CMakeLists.txt;0;")
add_test(TUV_5_4_COMPARE "python3" "test/json_yaml_compare.py" "example_tuv_5_4" "example_tuv_5_4_yaml")
set_tests_properties(TUV_5_4_COMPARE PROPERTIES  WORKING_DIRECTORY "/home/runner/work/tuv-x/tuv-x/build_test2" _BACKTRACE_TRIPLES "/home/runner/work/tuv-x/tuv-x/test/CMakeLists.txt;62;add_test;/home/runner/work/tuv-x/tuv-x/test/CMakeLists.txt;0;")
add_test(TS1_TSMLT_YAML "/home/runner/work/tuv-x/tuv-x/build_test2/tuv-x" "../examples/ts1_tsmlt.yml")
set_tests_properties(TS1_TSMLT_YAML PROPERTIES  WORKING_DIRECTORY "/home/runner/work/tuv-x/tuv-x/build_test2/example_ts1_tsmlt_yaml" _BACKTRACE_TRIPLES "/home/runner/work/tuv-x/tuv-x/test/CMakeLists.txt;68;add_test;/home/runner/work/tuv-x/tuv-x/test/CMakeLists.txt;0;")
add_test(TS1_TSMLT_COMPARE "python3" "test/json_yaml_compare.py" "example_ts1_tsmlt" "example_ts1_tsmlt_yaml")
set_tests_properties(TS1_TSMLT_COMPARE PROPERTIES  WORKING_DIRECTORY "/home/runner/work/tuv-x/tuv-x/build_test2" _BACKTRACE_TRIPLES "/home/runner/work/tuv-x/tuv-x/test/CMakeLists.txt;70;add_test;/home/runner/work/tuv-x/tuv-x/test/CMakeLists.txt;0;")
subdirs("unit")
