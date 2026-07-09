################################################################################
# Regression tests: quantum-yield algorithms

include(test_util)

create_regression_cxx_test(
    NAME    regression_quantum_yields
    SOURCES quantum_yields.cpp
)

################################################################################
