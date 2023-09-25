! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0

module tuvx_grid_equal_delta
! one dimension, equally spaced grid type. See
! :ref:`configuration-grids` for more information.

  use musica_constants,                only : dk => musica_dk
  use tuvx_grid,                       only : grid_t

  implicit none

  public :: grid_equal_delta_t

  type, extends(grid_t) :: grid_equal_delta_t
  contains
  end type grid_equal_delta_t

  !> Constructor
  interface grid_equal_delta_t
    module procedure constructor
  end interface grid_equal_delta_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function constructor( config ) result ( this )
    ! Initialize grid

    use musica_assert,                 only : assert_msg
    use musica_config,                 only : config_t
    use musica_string,                 only : string_t

    type(config_t), intent(inout) :: config ! The grid config. See :ref:`configuration-grids` for more details

    ! local variables
    integer  :: n
    real(dk) :: Lower_val, Upper_val, Delta_val
    character(len=*), parameter :: Iam = 'EqualDelta grid initialize: '
    type(grid_equal_delta_t), pointer  :: this
    type(string_t) :: required_keys(6), optional_keys(0)

    required_keys(1) = "type"
    required_keys(2) = "units"
    required_keys(3) = "begins at"
    required_keys(4) = "ends at"
    required_keys(5) = "cell delta"
    required_keys(6) = "name"

    call assert_msg( 127910504,                                               &
                     config%validate( required_keys, optional_keys ),         &
                     "Bad configuration data format for "//                   &
                     "equal interval grid." )

    allocate ( this )

    call config%get( 'begins at', Lower_val, Iam )
    call config%get( 'ends at', Upper_val, Iam )
    call config%get( 'cell delta', Delta_val, Iam )
    call config%get( 'name', this%handle_, Iam )
    call config%get( 'units', this%units_, Iam )

    this%ncells_ = ( Upper_val - Lower_val ) / Delta_val
    if( mod( ( Upper_val - Lower_val ), Delta_val ) /= 0._dk ) then
      this%ncells_ = this%ncells_ + 1
    endif
    allocate( this%mid_( this%ncells_ ) )
    allocate( this%delta_( this%ncells_ ) )
    allocate( this%edge_( this%ncells_ + 1 ) )
    do n = 1, this%ncells_ + 1
      this%edge_( n ) =                                                       &
        min( real( ( n - 1 ), kind = dk ) * Delta_val + Lower_val, Upper_val )
    enddo
    this%mid_(:) = .5_dk * &
      ( this%edge_( 1 : this%ncells_ ) + this%edge_( 2 : this%ncells_ + 1 ) )
    this%delta_(:) =                                                          &
      this%edge_( 2 : this%ncells_ + 1 ) - this%edge_( 1 : this%ncells_ )

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module tuvx_grid_equal_delta
