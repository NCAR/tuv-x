################################################################################
# Test utilities

include(test_util)

################################################################################
# Phototolysis rate regression tests

add_regression_test(
  regression_photolysis_rates
  test/regression/photolysis_rates/photolysis_rates.sh
  test/regression/photolysis_rates/photolysis_rates_memcheck.sh
)

################################################################################
