! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0

module tuvx_profile_earth_sun_distance
  ! Earth-Sun distance profile type
  ! Calculates the earth sun distance using the 
  ! :f:func:`tuvx_profile_utils/earth_sun_distance`

  use musica_constants,                only : dk => musica_dk
  use tuvx_profile,                    only : profile_t
  use tuvx_profile_utils,              only : julian_day_of_year,             &
                                              earth_sun_distance

  implicit none

  private
  public :: profile_earth_sun_distance_t

  type, extends(profile_t) :: profile_earth_sun_distance_t
  contains
    final     :: finalize
  end type profile_earth_sun_distance_t

  !> Constructor
  interface profile_earth_sun_distance_t
    module procedure constructor
  end interface profile_earth_sun_distance_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function constructor( config, grid_warehouse ) result ( this )
    ! Initialize distance between sun, earth in AU

    use musica_assert,                 only : assert_msg
    use musica_config,                 only : config_t
    use musica_string,                 only : string_t
    use tuvx_grid,                     only : grid_t
    use tuvx_grid_warehouse,           only : grid_warehouse_t

    ! Arguments
    type(profile_earth_sun_distance_t), pointer   :: this ! This f:type:`~tuvx_profile_earth_sun_distance/profile_earth_sun_distance_t`
    type(config_t), intent(inout)         :: config ! A profile config
    type(grid_warehouse_t), intent(inout) :: grid_warehouse ! A :f:type:`~tuvx_grid_warehouse/grid_warehouse_t`

    ! Local variables
    real(dk), parameter :: NINETY  = 90._dk
    integer :: tNdx
    integer :: Year, Month, Day
    integer :: Jday
    real(dk)    :: tmzone, ut, soldst
    character(len=*), parameter :: Iam = 'earth sun distance initialize: '
    class(grid_t), pointer :: timeGrid
    type(string_t) :: required_keys(6), optional_keys(1)

    required_keys(1) = "type"
    required_keys(2) = "units"
    required_keys(3) = "year"
    required_keys(4) = "month"
    required_keys(5) = "day"
    required_keys(6) = "name"
    optional_keys(1) = "time zone"

    call assert_msg( 236065474,                                               &
                     config%validate( required_keys, optional_keys ),         &
                     "Bad configuration data format for "//                   &
                     "Earth-Sun distance profile." )

    allocate ( this )

    timeGrid => grid_warehouse%get_grid( "time", "hours" )
    this%ncells_ = timeGrid%ncells_


    allocate( this%edge_val_(0) )

    !> Map solar zenith angle as function of time
    call config%get( 'name', this%handle_, Iam )
    call config%get( 'units', this%units_, Iam )
    call config%get( 'year', Year, Iam )
    call config%get( 'month', Month, Iam )
    call config%get( 'day', Day, Iam )
    call config%get( 'time zone', tmzone, Iam, default = 0.0_dk )

    Jday = julian_day_of_year( Year, Month, Day )

    do tNdx = 1, this%ncells_ + 1
      ut = timeGrid%edge_( tNdx ) - tmzone
      soldst = earth_sun_distance( Year, Jday, ut )
      this%edge_val_  = [ this%edge_val_, soldst ]
    enddo

    this%mid_val_ = .5_dk * ( this%edge_val_( 1 : this%ncells_ ) +            &
      this%edge_val_( 2 : this%ncells_ + 1 ) )

    this%delta_val_ = ( this%edge_val_( 2 : this%ncells_ + 1 ) -              &
      this%edge_val_( 1 : this%ncells_ ) )

    deallocate( timeGrid )

  end function constructor


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine finalize( this )
    ! Cleanup the memory used by this object

    type(profile_earth_sun_distance_t), intent(inout) :: this ! This f:type:`~tuvx_profile_earth_sun_distance/profile_earth_sun_distance_t`

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

end module tuvx_profile_earth_sun_distance
