! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: GPL-2.0-or-later
!
!vertical profile type
module micm_vert_Profile

  use musica_constants, only : musica_dk, musica_ik
  use musica_string,    only : string_t

  implicit none

  private
  public :: abs_vert_Profile_t, abs_vert_Profile_ptr

  type, abstract ::  abs_vert_Profile_t
    !> grid handle
    type(string_t) :: handle_
    !> number of wavelength grid cells
    integer(musica_ik)           :: ncells_
    !> scale heigth
    real(musica_dk)              :: hscale_
    !> cell centers
    real(musica_dk), allocatable :: mid_val_(:)
    !> cell edges
    real(musica_dk), allocatable :: edge_val_(:)
    !> cell deltas
    real(musica_dk), allocatable :: delta_val_(:)
    !> layer densities
    real(musica_dk), allocatable :: layer_dens_(:)
    !> layer densities including "exo" model layer
    real(musica_dk), allocatable :: exo_layer_dens_(:)
    !> overhead column burden
    real(musica_dk), allocatable :: burden_dens_(:)
  contains
    !> Initialize grid
    procedure(initial), deferred :: initialize
  end type abs_vert_Profile_t

  !> Pointer type for building sets of spectral wght objects
  type :: abs_vert_Profile_ptr
    class(abs_vert_Profile_t), pointer :: ptr_ => null( )
  end type abs_vert_Profile_ptr

interface

    !> Initialize grid
    subroutine initial( this, profile_config, gridWareHouse )
      
      use musica_config, only : config_t
      use musica_constants, only : musica_dk
      use micm_grid_warehouse,  only : grid_warehouse_t

      import abs_vert_Profile_t
      class(abs_vert_Profile_t), intent(inout) :: this
      type(config_t), intent(inout)            :: profile_config
      type(grid_warehouse_t), intent(inout)    :: gridWareHouse
    end subroutine initial

end interface

end module micm_vert_Profile
