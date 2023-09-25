! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0

module tuvx_profile_factory
  ! Build :f:type:`~tuvx_profile/profile_t` objects

  use tuvx_profile,                    only : profile_t
  use tuvx_profile_air,                only : profile_air_t
  use tuvx_profile_from_config,        only : profile_from_config_t
  use tuvx_profile_earth_sun_distance, only : profile_earth_sun_distance_t
  use tuvx_profile_extraterrestrial_flux,                                     &
    only : profile_extraterrestrial_flux_t
  use tuvx_profile_from_csv_file,      only : profile_from_csv_file_t
  use tuvx_profile_from_host,          only : profile_from_host_t
  use tuvx_profile_o2,                 only : profile_o2_t
  use tuvx_profile_o3,                 only : profile_o3_t
  use tuvx_profile_solar_zenith_angle, only : profile_solar_zenith_angle_t

  implicit none

  private
  public :: profile_builder, profile_type_name, profile_allocate

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function profile_builder( config, grid_warehouse ) result( new_profile_t )
    ! Build an instance of a :f:type:`~tuvx_profile/profile_t`

    use musica_assert,                   only : die_msg
    use musica_config,                   only : config_t
    use musica_string,                   only : string_t
    use tuvx_grid_warehouse,             only : grid_warehouse_t

    type(config_t),         intent(inout) :: config ! Grid configuration data
    type(grid_warehouse_t), intent(inout) :: grid_warehouse ! A :f:type:`~tuvx_grid_warehouse/grid_warehouse_t`

    class(profile_t), pointer :: new_profile_t ! New :f:type:`~tuvx_profile/profile_t` object

    ! Local variables
    character(len=*), parameter :: Iam = 'profile builder: '
    type(string_t) :: profile_type

    new_profile_t => null()
    call config%get( 'type', profile_type, Iam )

    select case( profile_type%to_char() )
      case( 'from csv file' )
        new_profile_t => profile_from_csv_file_t( config, grid_warehouse )
      case( 'extraterrestrial flux' )
        new_profile_t => profile_extraterrestrial_flux_t( config, grid_warehouse )
      case( 'from config file' )
        new_profile_t => profile_from_config_t( config, grid_warehouse )
      case( 'air' )
        new_profile_t => profile_air_t( config, grid_warehouse )
      case( 'O2' )
        new_profile_t => profile_o2_t( config, grid_warehouse )
      case( 'O3' )
        new_profile_t => profile_o3_t( config, grid_warehouse )
      case( 'solar zenith angle' )
        new_profile_t => profile_solar_zenith_angle_t( config, grid_warehouse )
      case( 'Earth-Sun distance' )
        new_profile_t => profile_earth_sun_distance_t( config, grid_warehouse )
      case default
        call die_msg( 884374015, "Invalid profile type: '" // &
          profile_type%to_char()//"'" )
    end select

  end function profile_builder

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  type(string_t) function profile_type_name( profile ) result( name )
    ! Returns the type of a profile as a string

    use musica_assert,                 only : die
    use musica_string,                 only : string_t

    class(profile_t), intent(in) :: profile ! profile to return type for

    select type( profile )
      type is( profile_from_csv_file_t )
        name = "profile_from_csv_file_t"
      type is( profile_extraterrestrial_flux_t )
        name = "profile_extraterrestrial_flux_t"
      type is( profile_from_config_t )
        name = "profile_from_config_t"
      type is( profile_from_host_t )
        name = "profile_from_host_t"
      type is( profile_air_t )
        name = "profile_air_t"
      type is( profile_o2_t )
        name = "profile_o2_t"
      type is( profile_o3_t )
        name = "profile_o3_t"
      type is( profile_solar_zenith_angle_t )
        name = "profile_solar_zenith_angle_t"
      type is( profile_earth_sun_distance_t )
        name = "profile_earth_sun_distance_t"
      class default
        call die( 616084530 )
      end select

  end function profile_type_name

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function profile_allocate( type_name ) result( profile )
    ! Allocates a profile pointer as a subclass by type

    use musica_assert,                 only : die_msg
    use musica_string,                 only : string_t

    type(string_t),   intent(in) :: type_name ! name of the type to allocate
    class(profile_t), pointer    :: profile   ! allocated profile

    select case( type_name%to_char( ) )
      case( 'profile_from_csv_file_t' )
        allocate( profile_from_csv_file_t :: profile )
      case( 'profile_extraterrestrial_flux_t' )
        allocate( profile_extraterrestrial_flux_t :: profile )
      case( 'profile_from_config_t' )
        allocate( profile_from_config_t :: profile )
      case( 'profile_from_host_t' )
        allocate( profile_from_host_t :: profile )
      case( 'profile_air_t' )
        allocate( profile_air_t :: profile )
      case( 'profile_o2_t' )
        allocate( profile_o2_t :: profile )
      case( 'profile_o3_t' )
        allocate( profile_o3_t :: profile )
      case( 'profile_solar_zenith_angle_t' )
        allocate( profile_solar_zenith_angle_t :: profile )
      case( 'profile_earth_sun_distance_t' )
        allocate( profile_earth_sun_distance_t :: profile )
      case default
        call die_msg( 706344703, "Invalid profile type: '"//type_name//"'" )
    end select

  end function profile_allocate

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module tuvx_profile_factory
