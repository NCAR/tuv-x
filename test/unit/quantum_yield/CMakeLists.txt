################################################################################
# Test utilities

include(test_util)

################################################################################
# Quantum yield tests

create_standard_test(NAME quantum_yield SOURCES base.F90 )
create_standard_test(NAME quantum_yield_no2_tint SOURCES no2_tint.F90 )
create_standard_test(NAME quantum_yield_tint SOURCES tint.F90 )

################################################################################
