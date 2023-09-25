################################################################################
# Radiative transfer source

add_subdirectory(radiators)
add_subdirectory(solvers)

target_sources(tuvx_object
  PRIVATE
    radiative_transfer.F90
    radiator.F90
    radiator_factory.F90
    radiator_warehouse.F90
    solver.F90
    solver_factory.F90
)

################################################################################
