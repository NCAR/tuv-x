################################################################################
# Test utilities

include(test_util)

################################################################################
# Phototolysis rate regression tests

# windows doesn't run .sh files
if(NOT WIN32)
  add_regression_test(
    regression_photolysis_rates
    test/regression/photolysis_rates/photolysis_rates.sh
    test/regression/photolysis_rates/photolysis_rates_memcheck.sh
  )
endif()

################################################################################
