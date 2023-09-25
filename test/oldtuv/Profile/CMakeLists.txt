################################################################################
# One dimensional grid source

set(PROFILE_SRC Profile/abstract.Profile.type.F90
                Profile/from_csv_file.type.F90
                Profile/from_config.type.F90
                Profile/srfAlbedo_from_config.type.F90
                Profile/air.from_csv_file.type.F90
                Profile/o2.from_csv_file.type.F90
                Profile/o3.from_csv_file.type.F90
                Profile/sza_from_time.type.F90
                Profile/earth_sun_distance.type.F90
                Profile/Profile_factory.F90
                Profile/Profile_warehouse.F90
            PARENT_SCOPE)

################################################################################
