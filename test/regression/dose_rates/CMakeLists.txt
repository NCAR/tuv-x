################################################################################
# Test utilities

include(test_util)

################################################################################
# Dose rate regression tests

# windows doesn't run .sh files
if(NOT WIN32)
  add_regression_test(
    regression_dose_rates
    test/regression/dose_rates/dose_rates.sh
    test/regression/dose_rates/dose_rates_memcheck.sh
  )
endif()

################################################################################
