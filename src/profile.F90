! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0

module tuvx_profile
  ! Profile type

  use musica_constants, only : musica_dk, musica_ik
  use musica_string,    only : string_t

  implicit none

  private
  public :: profile_t, profile_ptr

  type profile_t
    type(string_t) :: handle_ ! grid handle
    type(string_t) :: units_ ! units
    integer(musica_ik)           :: ncells_ ! number of grid cells
    real(musica_dk)              :: hscale_ ! scale height
    real(musica_dk), allocatable :: mid_val_(:) ! cell centers
    real(musica_dk), allocatable :: edge_val_(:) ! cell edges
    real(musica_dk), allocatable :: delta_val_(:) ! cell deltas
    real(musica_dk), allocatable :: layer_dens_(:) ! layer densities
    real(musica_dk), allocatable :: exo_layer_dens_(:) ! layer densities including "exo" model layer
    real(musica_dk), allocatable :: burden_dens_(:) ! overhead column burden
    logical                      :: enable_diagnostics ! determins if diagnostic output is written or not
  contains
    ! returns the units of a profile
    procedure :: units
    ! returns the number of grid cells for the profile
    procedure :: size => number_of_cells
    ! returns the number of bytes needed to pack the profile onto a buffer
    procedure :: pack_size
    ! packs the profile onto a character buffer
    procedure :: mpi_pack
    ! unpacks a profile from a character buffer into the object
    procedure :: mpi_unpack
    ! outputs profile data to a specified IO unit
    procedure :: output
  end type profile_t

  type profile_ptr
    ! Pointer type for building sets of profile objects
    class(profile_t), pointer :: val_ => null( )
  end type profile_ptr

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  type(string_t) function units( this )
    ! Returns the units for the profile

    class(profile_t), intent(in) :: this

    units = this%units_

  end function units

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function number_of_cells( this )
    ! Returns the number of grid cells for the profile

    class(profile_t), intent(in) :: this

    number_of_cells = this%ncells_

  end function number_of_cells

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function pack_size( this, comm )
    ! Returns the number of bytes required to pack the profile onto a buffer

    use musica_mpi,                    only : musica_mpi_pack_size

    class(profile_t),  intent(in) :: this ! profile to be packed
    integer,           intent(in) :: comm ! MPI communicator

#ifdef MUSICA_USE_MPI
    pack_size = this%handle_%pack_size( comm )                                &
                + this%units_%pack_size( comm )                               &
                + musica_mpi_pack_size( this%ncells_,         comm )          &
                + musica_mpi_pack_size( this%hscale_,         comm )          &
                + musica_mpi_pack_size( this%mid_val_,        comm )          &
                + musica_mpi_pack_size( this%edge_val_,       comm )          &
                + musica_mpi_pack_size( this%delta_val_,      comm )          &
                + musica_mpi_pack_size( this%layer_dens_,     comm )          &
                + musica_mpi_pack_size( this%exo_layer_dens_, comm )          &
                + musica_mpi_pack_size( this%burden_dens_,    comm )
#else
    pack_size = 0
#endif

  end function pack_size

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine mpi_pack( this, buffer, position, comm )
    ! Packs the profile onto a character buffer

    use musica_assert,                 only : assert
    use musica_mpi,                    only : musica_mpi_pack

    class(profile_t),  intent(in)    :: this      ! profile to be packed
    character,         intent(inout) :: buffer(:) ! memory buffer
    integer,           intent(inout) :: position  ! current buffer position
    integer,           intent(in)    :: comm      ! MPI communicator

#ifdef MUSICA_USE_MPI
    integer :: prev_pos

    prev_pos = position
    call this%handle_%mpi_pack( buffer, position, comm )
    call this%units_%mpi_pack(  buffer, position, comm )
    call musica_mpi_pack( buffer, position, this%ncells_,         comm )
    call musica_mpi_pack( buffer, position, this%hscale_,         comm )
    call musica_mpi_pack( buffer, position, this%mid_val_,        comm )
    call musica_mpi_pack( buffer, position, this%edge_val_,       comm )
    call musica_mpi_pack( buffer, position, this%delta_val_,      comm )
    call musica_mpi_pack( buffer, position, this%layer_dens_,     comm )
    call musica_mpi_pack( buffer, position, this%exo_layer_dens_, comm )
    call musica_mpi_pack( buffer, position, this%burden_dens_,    comm )
    call assert( 947035840, position - prev_pos <= this%pack_size( comm ) )
#endif

  end subroutine mpi_pack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine mpi_unpack( this, buffer, position, comm )
    ! Unpacks the profile from a character buffer

    use musica_assert,                 only : assert
    use musica_mpi,                    only : musica_mpi_unpack

    class(profile_t),  intent(out)   :: this      ! profile to be unpacked
    character,         intent(inout) :: buffer(:) ! memory buffer
    integer,           intent(inout) :: position  ! current buffer position
    integer,           intent(in)    :: comm      ! MPI communicator

#ifdef MUSICA_USE_MPI
    integer :: prev_pos

    prev_pos = position
    call this%handle_%mpi_unpack( buffer, position, comm )
    call this%units_%mpi_unpack(  buffer, position, comm )
    call musica_mpi_unpack( buffer, position, this%ncells_,         comm )
    call musica_mpi_unpack( buffer, position, this%hscale_,         comm )
    call musica_mpi_unpack( buffer, position, this%mid_val_,        comm )
    call musica_mpi_unpack( buffer, position, this%edge_val_,       comm )
    call musica_mpi_unpack( buffer, position, this%delta_val_,      comm )
    call musica_mpi_unpack( buffer, position, this%layer_dens_,     comm )
    call musica_mpi_unpack( buffer, position, this%exo_layer_dens_, comm )
    call musica_mpi_unpack( buffer, position, this%burden_dens_,    comm )
    call assert( 137296014, position - prev_pos <= this%pack_size( comm ) )
#endif

  end subroutine mpi_unpack

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine output( this, grid, io_unit )
    ! Outputs profile data to a specified IO unit (or stdout by default)

    use musica_assert,                 only : assert_msg
    use musica_string,                 only : to_char
    use tuvx_grid,                     only : grid_t

    class(profile_t),        intent(in) :: this
    class(grid_t), optional, intent(in) :: grid ! grid to plot values against
    integer,       optional, intent(in) :: io_unit

    integer :: io, i_cell
    type(string_t) :: axis_label, exo_val, layer_val, burden_val
    real(kind=musica_dk) :: axis_val

    io = 6
    if( present( io_unit ) ) io = io_unit
    layer_val = "---"
    exo_val = "---"
    burden_val = "---"
    write(io,*) "# Profile: "//this%handle_%val_//" ("//this%units_%val_//")"
    write(io,*) "# number of cells:", this%ncells_
    if( present( grid ) ) then
      axis_label%val_ = grid%handle_%val_//" ("//grid%units_%val_//")"
      call assert_msg( 207788261, grid%ncells_ .eq. this%ncells_,             &
                       "Grid cell mismatch printing profile. Expected "//     &
                       trim( to_char( this%ncells_ ) )//" but got "//         &
                       trim( to_char( this%ncells_ ) ) )
    else
      axis_label%val_ = "index"
    end if
    write(io,*) "# "//axis_label%val_//", "//                                 &
          "mid-point, delta, layer density, exo layer density, burden density"
    do i_cell = 1, this%ncells_
      axis_val = real( i_cell, kind=musica_dk )
      if( present( grid ) ) axis_val = grid%mid_( i_cell )
      if( allocated( this%layer_dens_ ) ) then
        layer_val%val_ = trim( to_char( this%layer_dens_( i_cell ) ) )
      end if
      if( allocated( this%exo_layer_dens_ ) ) then
        exo_val%val_ = trim( to_char( this%exo_layer_dens_( i_cell ) ) )
      end if
      if( allocated( this%burden_dens_ ) ) then
        burden_val%val_ = trim( to_char( this%burden_dens_( i_cell ) ) )
      end if
      write(io,*) axis_val, ",", this%mid_val_( i_cell ), ",",                &
                  this%delta_val_( i_cell ), ",", layer_val%val_, ",",        &
                  exo_val%val_, ",", burden_val%val_
    end do
    if( allocated( this%exo_layer_dens_ ) ) then
      write(io,*) "---, ---, ---, ---,",                                      &
                  this%exo_layer_dens_( this%ncells_ + 1 ), ", ---"
    end if
    write(*,*) "# "//axis_label%val_//", edge"
    do i_cell = 1, this%ncells_ + 1
      axis_val = real( i_cell, kind=musica_dk )
      if( present( grid ) ) axis_val = grid%edge_( i_cell )
      write(io,*) axis_val, ",", this%edge_val_( i_cell )
    end do

  end subroutine output

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end module tuvx_profile
