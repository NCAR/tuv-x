! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
  
module tuvx_profile_solar_zenith_angle
  ! Solar zenith angle from time.
  ! Calculates the solar zenith angle using the 
  ! :f:func:`tuvx_profile_utils/solar_zenith_angle`

  use musica_constants,   only : &
    dk => musica_dk, ik => musica_ik, lk => musica_lk
  use tuvx_profile,       only : profile_t
  use musica_assert,      only : die_msg
  use tuvx_profile_utils, only : julian_day_of_year, solar_zenith_angle

  implicit none

  public :: profile_solar_zenith_angle_t

  type, extends(profile_t) :: profile_solar_zenith_angle_t
  contains
    final     :: finalize
  end type profile_solar_zenith_angle_t

  interface profile_solar_zenith_angle_t
    module procedure constructor
  end interface profile_solar_zenith_angle_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function constructor( config, grid_warehouse ) result( this )
    ! Initialize solar zenith angle profile from time

    use musica_assert,                 only : assert_msg
    use musica_config,                 only : config_t
    use musica_string,                 only : string_t
    use tuvx_grid,                     only : grid_t
    use tuvx_grid_warehouse,           only : grid_warehouse_t

    type(profile_solar_zenith_angle_t), pointer   :: this ! This f:type:`~tuvx_profile_solar_zenith_angle/profile_solar_zenith_angle_t`
    type(config_t), intent(inout)         :: config ! A profile config
    type(grid_warehouse_t), intent(inout) :: grid_warehouse ! A :f:type:`~tuvx_grid_warehouse/grid_warehouse_t`

    ! Local variables
    real(dk), parameter :: NINETY  = 90._dk
    integer(ik) :: tNdx
    integer(ik) :: Year, Month, Day
    integer(ik) :: Jday
    real(dk)    :: tmzone, ut, solarElevation
    real(dk)    :: Lon, Lat
    character(len=*), parameter :: Iam = 'sza from time initialize: '
    class(grid_t), pointer :: timeGrid
    type(string_t) :: required_keys(8), optional_keys(1)

    required_keys(1) = "type"
    required_keys(2) = "units"
    required_keys(3) = "year"
    required_keys(4) = "month"
    required_keys(5) = "day"
    required_keys(6) = "longitude"
    required_keys(7) = "latitude"
    required_keys(8) = "name"
    optional_keys(1) = "time zone"

    call assert_msg( 482373225,                                               &
                     config%validate( required_keys, optional_keys ),         &
                     "Bad configuration data format for "//                   &
                     "solar zenith angle profile." )
    allocate( this )

    timeGrid => grid_warehouse%get_grid( "time", "hours" )
    this%ncells_ = timeGrid%ncells_


    allocate( this%edge_val_(0) )

    ! Map solar zenith angle as function of time
    call config%get( 'name', this%handle_, Iam )
    call config%get( 'units', this%units_, Iam )
    call config%get( 'year', Year, Iam )
    call config%get( 'month', Month, Iam )
    call config%get( 'day', Day, Iam )
    call config%get( 'time zone', tmzone, Iam, default=0.0_dk )
    call config%get( 'longitude', Lon, Iam )
    call config%get( 'latitude', Lat, Iam )

    Jday = julian_day_of_year(Year, Month, Day )

    do tNdx = 1_ik,this%ncells_+1_ik
      ut = timeGrid%edge_(tNdx) + tmzone
      solarElevation = solar_zenith_angle(Year, Jday, ut, Lat, Lon )
      this%edge_val_ = [this%edge_val_,NINETY - solarElevation]
    enddo

    this%mid_val_ = .5_dk * ( &
      this%edge_val_(1_ik:this%ncells_) + &
      this%edge_val_(2_ik:this%ncells_+1_ik) &
    )

    this%delta_val_ = this%edge_val_(2_ik:this%ncells_+1_ik) - &
      this%edge_val_(1_ik:this%ncells_)

    deallocate( timeGrid )

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine finalize( this )
    ! Cleanup the memory used by this object

    type(profile_solar_zenith_angle_t), intent(inout) :: this ! This f:type:`~tuvx_profile_solar_zenith_angle/profile_solar_zenith_angle_t`

    if( allocated( this%edge_val_ ) ) then
      deallocate( this%edge_val_ )
    endif
    if( allocated( this%mid_val_ ) ) then
      deallocate( this%mid_val_ )
    endif
    if( allocated( this%delta_val_ ) ) then
      deallocate( this%delta_val_ )
    endif

  end subroutine finalize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module tuvx_profile_solar_zenith_angle
