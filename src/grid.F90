! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
module tuvx_grid
! A one dimensional grid type.
!

  use musica_constants, only : musica_dk, musica_ik
  use musica_string,    only : string_t

  implicit none

  private
  public :: grid_t, grid_ptr

  type ::  grid_t
    type(string_t) :: handle_ ! grid handle
    type(string_t) :: units_ ! units
    integer(musica_ik) :: ncells_ ! number of wavelength grid cells
    real(musica_dk), allocatable :: mid_(:) ! cell centers
    real(musica_dk), allocatable :: edge_(:) ! cell edges
    real(musica_dk), allocatable :: delta_(:) ! cell deltas
  contains
    ! Returns the number of grid cells
    procedure :: size => number_of_cells
    ! Returns the units for the grid
    procedure :: units
    ! Returns the number of bytes needed to pack the grid onto a buffer
    procedure :: pack_size
    ! Packs the grid onto a character buffer
    procedure :: mpi_pack
    ! Unpacks a grid from a character buffer into the object
    procedure :: mpi_unpack
    ! Output the grid state
    procedure :: output
  end type grid_t

  !> Pointer type for building sets of spectral wght objects
  type :: grid_ptr
    class(grid_t), pointer :: val_ => null( )
  end type grid_ptr

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function number_of_cells( this )
    ! Returns the number of grid cells

    class(grid_t), intent(in) :: this ! grid

    number_of_cells = this%ncells_

  end function number_of_cells

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  type(string_t) function units( this )
  ! Returns the units for the grid

    class(grid_t), intent(in) :: this ! A :f:type:`~tuvx_grid/grid_t`

    units = this%units_

  end function units

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function pack_size( this, comm )
    ! Returns the number of bytes needed to pack the grid onto a buffer

    use musica_mpi,                    only : musica_mpi_pack_size

    class(grid_t),     intent(in) :: this ! grid to pack
    integer,           intent(in) :: comm ! MPI communicator

#ifdef MUSICA_USE_MPI
    pack_size = this%handle_%pack_size( comm ) +                              &
                this%units_%pack_size( comm ) +                               &
                musica_mpi_pack_size( this%ncells_, comm ) +                  &
                musica_mpi_pack_size( this%mid_,    comm ) +                  &
                musica_mpi_pack_size( this%edge_,   comm ) +                  &
                musica_mpi_pack_size( this%delta_,  comm )
#else
    pack_size = 0
#endif

  end function pack_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine mpi_pack( this, buffer, position, comm )
    ! Packs the grid onto a character buffer

    use musica_assert,                 only : assert
    use musica_mpi,                    only : musica_mpi_pack

    class(grid_t),     intent(in)    :: this      ! grid to pack
    character,         intent(inout) :: buffer(:) ! memory buffer
    integer,           intent(inout) :: position  ! current buffer
    integer,           intent(in)    :: comm      ! MPI communicator

#ifdef MUSICA_USE_MPI
    integer :: prev_pos

    prev_pos = position
    call this%handle_%mpi_pack( buffer, position, comm )
    call this%units_%mpi_pack(  buffer, position, comm )
    call musica_mpi_pack( buffer, position, this%ncells_, comm )
    call musica_mpi_pack( buffer, position, this%mid_,     comm )
    call musica_mpi_pack( buffer, position, this%edge_,    comm )
    call musica_mpi_pack( buffer, position, this%delta_,   comm )
    call assert( 887266509, position - prev_pos <= this%pack_size( comm ) )
#endif

  end subroutine mpi_pack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine mpi_unpack( this, buffer, position, comm )
    ! Unpacks the grid from a character buffer into the object

    use musica_assert,                 only : assert
    use musica_mpi,                    only : musica_mpi_unpack

    class(grid_t),     intent(out)   :: this      ! grid to be unpacked
    character,         intent(inout) :: buffer(:) ! memory buffer
    integer,           intent(inout) :: position  ! current buffer position
    integer,           intent(in)    :: comm      ! MPI communicator

#ifdef MUSICA_USE_MPI
    integer :: prev_pos

    prev_pos = position
    call this%handle_%mpi_unpack( buffer, position, comm )
    call this%units_%mpi_unpack(  buffer, position, comm )
    call musica_mpi_unpack( buffer, position, this%ncells_, comm )
    call musica_mpi_unpack( buffer, position, this%mid_,     comm )
    call musica_mpi_unpack( buffer, position, this%edge_,    comm )
    call musica_mpi_unpack( buffer, position, this%delta_,   comm )
    call assert( 805916539, position - prev_pos <= this%pack_size( comm ) )
#endif

  end subroutine mpi_unpack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine output( this, io_unit )
    ! Prints the grid state to the specified unit (stdout by default)

    class(grid_t),     intent(in) :: this
    integer, optional, intent(in) :: io_unit

    integer :: i_cell, io

    io = 6
    if( present( io_unit ) ) io = io_unit
    write(io,*) "# Grid: "//this%handle_%val_//" ("//this%units_%val_//")"
    write(io,*) "# Number of cells", this%ncells_
    write(io,*) "# index, mid-point, edge, width"
    do i_cell = 1, this%ncells_
      write(io,*) i_cell, ",", this%mid_( i_cell ), ",",                      &
                  this%edge_( i_cell ), ",", this%delta_( i_cell )
    end do
    write(io,*) this%ncells_ + 1, ", ---,", this%edge_( this%ncells_ + 1 ),   &
                ", ---"

  end subroutine output

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module tuvx_grid
