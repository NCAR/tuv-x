################################################################################
# utilities for cross section parameterizations

target_sources(tuvx_object
  PRIVATE
    temperature_parameterization.F90
    temperature_parameterization_burkholder.F90
    temperature_parameterization_harwood.F90
    temperature_parameterization_taylor_series.F90
    temperature_range.F90  
)

################################################################################