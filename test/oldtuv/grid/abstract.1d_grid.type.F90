! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: GPL-2.0-or-later
!
!one dimension grid type
module micm_1d_grid

  use musica_constants, only : musica_dk, musica_ik
  use musica_string,    only : string_t

  implicit none

  private
  public :: base_grid_t, base_grid_ptr

  type, abstract ::  base_grid_t
    !> grid handle
    type(string_t) :: handle_
    !> number of wavelength grid cells
    integer(musica_ik) :: ncells_
    !> cell centers
    real(musica_dk), allocatable :: mid_(:)
    !> cell edges
    real(musica_dk), allocatable :: edge_(:)
    !> cell deltas
    real(musica_dk), allocatable :: delta_(:)
  contains
    !> Initialize grid
    procedure(initial), deferred :: initialize
  end type base_grid_t

  !> Pointer type for building sets of spectral wght objects
  type :: base_grid_ptr
    class(base_grid_t), pointer :: ptr_ => null( )
  end type base_grid_ptr

interface

    !> Initialize grid
    subroutine initial( this, grid_config )
      
      use musica_config, only : config_t

      import base_grid_t
      class(base_grid_t), intent(inout) :: this
      type(config_t), intent(inout)       :: grid_config
    end subroutine initial

end interface

end module micm_1d_grid
