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
    integer(musica_ik) :: ncells_
    !> cell centers
    real(musica_dk), allocatable :: mid_val_(:)
    !> cell edges
    real(musica_dk), allocatable :: edge_val_(:)
    !> cell deltas
    real(musica_dk), allocatable :: delta_val_(:)
  contains
    !> Initialize grid
    procedure(initial), deferred :: initialize
    procedure                    :: inter1
  end type abs_vert_Profile_t

  !> Pointer type for building sets of spectral wght objects
  type :: abs_vert_Profile_ptr
    class(abs_vert_Profile_t), pointer :: ptr_ => null( )
  end type abs_vert_Profile_ptr

interface

    !> Initialize grid
    subroutine initial( this, profile_config, zGrid )
      
      use musica_config, only : config_t
      use musica_constants, only : musica_dk

      import abs_vert_Profile_t
      class(abs_vert_Profile_t), intent(inout) :: this
      type(config_t), intent(inout)       :: profile_config
      real(musica_dk), intent(in)         :: zGrid(:)
    end subroutine initial

end interface

    contains

    FUNCTION inter1(this, xtarget, xsrc,ysrc) result( ytarget )

    IMPLICIT NONE

    !> Arguments
    REAL(musica_dk), intent(in) :: xtarget(:)
    REAL(musica_dk), intent(in) :: xsrc(:), ysrc(:)

    REAL(musica_dk) :: ytarget(size(xtarget))
    class(abs_vert_Profile_t)    :: this

    !> Local variables
    INTEGER(musica_ik) :: n
    INTEGER(musica_ik) :: jsave, i, j
    REAL(musica_dk)    :: slope

    n = size(xsrc)
    jsave = 1
    ytarget = 0.
    DO i = 1,size(xtarget)
      do j = jsave,n-1
        IF( (xsrc(j) <= xtarget(i)) .and. (xtarget(i) <= xsrc(j+1)) ) THEN
          slope = (ysrc(j+1)-ysrc(j)) / (xsrc(j+1)-xsrc(j))
          ytarget(i) = ysrc(j) + slope * (xtarget(i) - xsrc(j))
          jsave = j
          exit
        ENDIF
      enddo
    ENDDO

    END FUNCTION inter1

end module micm_vert_Profile
