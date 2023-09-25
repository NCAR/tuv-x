! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
module tuvx_grid_from_host
! 1d grid whose data will be provided by the host application at runtime
! See :ref:`configuration-grids` for more information.

  use musica_constants,                only : dk => musica_dk
  use tuvx_grid,                       only : grid_t

  implicit none

  public :: grid_from_host_t, grid_updater_t

  type, extends(grid_t) :: grid_from_host_t
    ! grid that can be updated from a host application
  contains
  end type grid_from_host_t

  ! grid constructor
  interface grid_from_host_t
    module procedure :: constructor_char
    module procedure :: constructor_string
  end interface grid_from_host_t

  type :: grid_updater_t
#ifndef MUSICA_IS_NAG_COMPILER
    private
#endif
    ! updater for `grid_from_host_t` grids
    class(grid_from_host_t), pointer :: grid_ => null( )
  contains
    procedure :: update
  end type grid_updater_t

  ! grid updater constructor
  interface grid_updater_t
    module procedure :: updater_constructor
  end interface grid_updater_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function constructor_char( name, units, number_of_cells ) result( this )
    ! Initialize the grid

    use musica_assert,                 only : assert_msg
    use musica_string,                 only : to_char

    character(len=*),       intent(in) :: name            ! name of the grid
    character(len=*),       intent(in) :: units           ! units for the grid
    integer,                intent(in) :: number_of_cells ! number of grid cells
    type(grid_from_host_t), pointer    :: this            ! constructed grid

    allocate( this )

    this%handle_ = name
    this%units_  = units
    this%ncells_ = number_of_cells

    call assert_msg( 968598222, this%ncells_ >= 0,                            &
                     "Invalid grid size for grid from host: "//               &
                     trim( to_char( number_of_cells ) ) )
    allocate( this%mid_(   this%ncells_     ) )
    allocate( this%edge_(  this%ncells_ + 1 ) )
    allocate( this%delta_( this%ncells_     ) )
    this%mid_(:)   = 0.0_dk
    this%edge_(:)  = 0.0_dk
    this%delta_(:) = 0.0_dk

  end function constructor_char

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function constructor_string( name, units, number_of_cells ) result( this )
    ! Initialize the grid

    use musica_string,                 only : string_t

    type(string_t),         intent(in) :: name            ! name of the grid
    type(string_t),         intent(in) :: units           ! units for the grid
    integer,                intent(in) :: number_of_cells ! number of grid cells
    type(grid_from_host_t), pointer    :: this            ! constructed grid

    this => constructor_char( name%to_char( ), units%to_char( ),              &
                              number_of_cells )

  end function constructor_string

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function updater_constructor( grid ) result( this )
    ! Constructs an updater for a `grid_from_host_t` grid

    class(grid_from_host_t), target, intent(inout) :: grid ! grid to be updated by host
    type(grid_updater_t)                           :: this ! new updater

    this%grid_ => grid

  end function updater_constructor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine update( this, edges, mid_points )
    ! Updates the target grid

    use musica_assert,                 only : assert_msg, die_msg
    use musica_string,                 only : to_char

    class(grid_updater_t),   intent(inout) :: this          ! grid updater
    real(kind=dk), optional, intent(in)    :: edges(:)      ! new edge values
    real(kind=dk), optional, intent(in)    :: mid_points(:) ! new mid-point values

    integer :: size_grid, size_host

    call assert_msg( 689055048, associated( this%grid_ ),                     &
                     "Cannot update an unspecified grid." )
    call assert_msg( 379490048, present( edges ),                             &
                     "Edges must be provided for grid update." )
    size_grid = size( this%grid_%edge_ )
    size_host = size( edges )
    if( size_grid .ne. size_host ) then
      call die_msg( 625263958,                                                &
                    "Size mismatch for grid edges for grid '"//               &
                    this%grid_%handle_//"'. Expected "//                      &
                    trim( to_char( size_grid ) )//", got "//                  &
                    trim( to_char( size_host ) ) )
    end if
    this%grid_%edge_(:) = edges(:)
    this%grid_%delta_(:) = this%grid_%edge_( 2 : this%grid_%ncells_ + 1 ) -   &
                           this%grid_%edge_( 1 : this%grid_%ncells_ )
    if( present( mid_points ) ) then
      size_grid = size( this%grid_%mid_ )
      size_host = size( mid_points )
      if( size_grid .ne. size_host ) then
        call die_msg( 742081595,                                              &
                      "Size mismatch for grid mid-points for grid '"//        &
                      this%grid_%handle_//"'. Expected "//                    &
                      trim( to_char( size_grid ) )//", got "//                &
                      trim( to_char( size_host ) ) )
      end if
      this%grid_%mid_(:) = mid_points(:)
    else
      this%grid_%mid_( 1 : this%grid_%ncells_ ) =                             &
        ( this%grid_%edge_( 1 : this%grid_%ncells_ ) +                        &
          this%grid_%edge_( 2 : this%grid_%ncells_ + 1 ) ) * 0.5_dk
    end if

  end subroutine update

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module tuvx_grid_from_host
