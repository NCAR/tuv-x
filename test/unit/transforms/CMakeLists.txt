################################################################################
# Transform system unit tests

include(test_util)

create_standard_cxx_test(NAME transforms_factories   SOURCES test_factories.cpp)
create_standard_cxx_test(NAME transforms_combinators SOURCES test_combinators.cpp)

################################################################################
