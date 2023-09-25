################################################################################
# Test utilities

include(test_util)

################################################################################
# Radiator regression tests

add_regression_test(
  regression_aerosols_radiator
  test/regression/radiators/radiation.aerosols.sh
  test/regression/radiators/radiation.aerosols.memcheck.sh
)

add_regression_test(
  regression_aerosols_radiator_4strm
  test/regression/radiators/radiation.aerosols.4strm.sh
  test/regression/radiators/radiation.aerosols.4strm.memcheck.sh
)

add_regression_test(
  regression_o2_radiator
  test/regression/radiators/radiation.o2.sh
  test/regression/radiators/radiation.o2.memcheck.sh
)

add_regression_test(
  regression_o2_radiator_4strm
  test/regression/radiators/radiation.o2.4strm.sh
  test/regression/radiators/radiation.o2.4strm.memcheck.sh
)

add_regression_test(
  regression_o3_radiator
  test/regression/radiators/radiation.o3.sh
  test/regression/radiators/radiation.o3.memcheck.sh
)

add_regression_test(
  regression_o3_radiator_4strm
  test/regression/radiators/radiation.o3.4strm.sh
  test/regression/radiators/radiation.o3.4strm.memcheck.sh
)

add_regression_test(
  regression_rayleigh_radiator
  test/regression/radiators/radiation.rayleigh.sh
  test/regression/radiators/radiation.rayleigh.memcheck.sh
)

add_regression_test(
  regression_rayleigh_radiator_4strm
  test/regression/radiators/radiation.rayleigh.4strm.sh
  test/regression/radiators/radiation.rayleigh.4strm.memcheck.sh
)

add_regression_test(
  regression_all_radiators
  test/regression/radiators/radiation.all.sh
  test/regression/radiators/radiation.all.memcheck.sh
)

add_regression_test(
  regression_all_radiators_4strm
  test/regression/radiators/radiation.all.4strm.sh
  test/regression/radiators/radiation.all.4strm.memcheck.sh
)

################################################################################
