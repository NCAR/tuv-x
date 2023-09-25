! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
module tuvx_grid_from_config
! 1d grid specified in json config file. See
! :ref:`configuration-grids` for more information.

  use musica_constants,                only : dk => musica_dk
  use tuvx_grid,                       only : grid_t

  implicit none

  public :: grid_from_config_t

  type, extends(grid_t) :: grid_from_config_t
  contains
  end type grid_from_config_t

  !> Constructor
  interface grid_from_config_t
    module procedure constructor
  end interface grid_from_config_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function constructor( config ) result ( this )
    ! Initialize grid

    use musica_assert,                 only : assert_msg
    use musica_config,                 only : config_t
    use musica_string,                 only : string_t

    type(config_t), intent(inout) :: config ! The grid config. See :ref:`configuration-grids` for more details

    ! Local variables
    character(len=*), parameter :: Iam = 'From config grid initialize: '
    type(grid_from_config_t), pointer  :: this
    type(string_t) :: required_keys(4), optional_keys(0)

    required_keys(1) = "type"
    required_keys(2) = "units"
    required_keys(3) = "values"
    required_keys(4) = "name"

    call assert_msg( 482373225,                                               &
                     config%validate( required_keys, optional_keys ),         &
                     "Bad configuration data format for "//                   &
                     "grid from configuration file." )
    allocate( this )

    call config%get( 'name', this%handle_, Iam )
    call config%get( 'units', this%units_, Iam )
    call config%get( "values", this%edge_, Iam )

    this%ncells_ = size( this%edge_ ) - 1
    this%mid_ = .5_dk *                                                       &
      ( this%edge_( 1 : this%ncells_ ) + this%edge_( 2 : this%ncells_ + 1 ) )
    this%delta_ =                                                             &
      this%edge_( 2 : this%ncells_ + 1 ) - this%edge_( 1 : this%ncells_ )

  end function constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module tuvx_grid_from_config
