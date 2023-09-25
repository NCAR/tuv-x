################################################################################
# Profile source

target_sources(tuvx_object
  PRIVATE
    air.F90
    earth_sun_distance.F90
    extraterrestrial_flux.F90
    from_config.F90
    from_csv_file.F90
    from_host.F90
    o2.F90
    o3.F90
    profile_utils.F90
    solar_zenith_angle.F90
)

################################################################################
