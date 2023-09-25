################################################################################
# Spectral weight source

target_sources(tuvx_object
  PRIVATE
    notch_filter.F90
    gaussian.F90
    eppley.F90
    par.F90
    exp_decay.F90
    scup_mice.F90
    standard_human_erythema.F90
    UV_Index.F90
    phytoplankton_boucher.F90
    plant_damage.F90
    plant_damage_flint_caldwell.F90
    plant_damage_flint_caldwell_ext.F90
)

################################################################################
