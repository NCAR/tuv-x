################################################################################
# Spectral weights source

set(SW_SRC spectral_wght/UV_Index.spectral_wght.type.F90
           spectral_wght/abstract.spectral_wght.type.F90
           spectral_wght/base.spectral_wght.type.F90
           spectral_wght/eppley_uv_photometer.spectral_wght.type.F90
           spectral_wght/exponential_decay.spectral_wght.type.F90
           spectral_wght/gaussian_305_nm_10_nm_FWHM.spectral_wght.type.F90
           spectral_wght/gaussian_320_nm_10_nm_FWHM.spectral_wght.type.F90
           spectral_wght/gaussian_340_nm_10_nm_FWHM.spectral_wght.type.F90
           spectral_wght/gaussian_380_nm_10_nm_FWHM.spectral_wght.type.F90
           spectral_wght/par_400-700nm.spectral_wght.type.F90
           spectral_wght/phytoplankton_boucher.spectral_wght.type.F90
           spectral_wght/plant_damage.spectral_wght.type.F90
           spectral_wght/plant_damage_flint_caldwell.spectral_wght.type.F90
           spectral_wght/plant_damage_flint_caldwell_ext.spectral_wght.type.F90
           spectral_wght/scup_mice.spectral_wght.type.F90
           spectral_wght/standard_human_erythema.spectral_wght.type.F90
           spectral_wght/uv-a_315_400_nm.spectral_wght.type.F90
           spectral_wght/uv-b_280_315_nm.spectral_wght.type.F90
           spectral_wght/uv-b_280_320_nm.spectral_wght.type.F90
           spectral_wght/visplus.spectral_wght.type.F90
           spectral_wght/spectral_wght_warehouse.F90
           spectral_wght/spectral_wght_factory.F90
           PARENT_SCOPE)

################################################################################
