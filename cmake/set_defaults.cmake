# Overwrite the init flags chosen by CMake
if(${CMAKE_Fortran_COMPILER_ID} MATCHES "GNU")
  set(CMAKE_Fortran_FLAGS_DEBUG_INIT "-g -O0 -fcheck=bounds,do,pointer -ffpe-trap=zero,overflow,invalid")
  set(CMAKE_Fortran_FLAGS_COVERAGE_INIT "g -O0 -fprofile-arcs -ftest-coverage -fprofile-abs-path -fcheck=bounds,do,pointer -ffpe-trap=zero,overflow,invalid")
elseif(${CMAKE_Fortran_COMPILER_ID} MATCHES "Intel")
  set(CMAKE_Fortran_FLAGS_DEBUG_INIT "-g -traceback -ip -O0 -fp-model precise -w -ftz -align all -fno-alias -convert big_endian")
endif()
