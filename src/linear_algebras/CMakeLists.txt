################################################################################
# Linear algebra source

target_sources(tuvx_object
  PRIVATE
    linpack.F90
)

if(TUVX_ENABLE_LAPACK)
  target_sources(tuvx_object PRIVATE lapack.F90)
endif()

################################################################################
