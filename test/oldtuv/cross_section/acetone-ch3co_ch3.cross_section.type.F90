! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: GPL-2.0-or-later
!
!> \file
!> This ch3coch3+hv->ch3co_ch3 cross_section module

!> The ch3coch3+hv->ch3co+ch3_cross_section type and related functions
module micm_ch3coch3_ch3co_ch3_cross_section_type

  use micm_base_cross_section_type,    only : base_cross_section_t
  use musica_constants,                only : musica_dk, musica_ik

  implicit none

  private
  public :: ch3coch3_ch3co_ch3_cross_section_t

  !> Calculator for acetone cross_section
  type, extends(base_cross_section_t) :: ch3coch3_ch3co_ch3_cross_section_t
  contains
    !> Initialize the cross section
    procedure :: calculate => run
  end type ch3coch3_ch3co_ch3_cross_section_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculate the photorate cross section for a given set of environmental conditions
  function run( this, environment ) result( cross_section )
!     qyacet - q.y. for acetone, based on Blitz et al. (2004)
! Compute acetone quantum yields according to the parameterization of:
! Blitz, M. A., D. E. Heard, M. J. Pilling, S. R. Arnold, and M. P. Chipperfield 
!       (2004), Pressure and temperature-dependent quantum yields for the 
!       photodissociation of acetone between 279 and 327.5 nm, Geophys. 
!       Res. Lett., 31, L06111, doi:10.1029/2003GL018793.

    use micm_environment,  only : environment_t
    use musica_assert,     only : die_msg

    !> base cross section
    class(ch3coch3_ch3co_ch3_cross_section_t), intent(in) :: this
    !> Environmental conditions
    class(environment_t), intent(in)             :: environment
    !> Calculated cross section
    real(kind=musica_dk)                         :: cross_section(size(this%mdl_lambda_center))

    real(musica_dk), parameter :: rONE  = 1.0_musica_dk
    character(len=*), parameter :: Iam = 'ch3coch3+hv->ch3co+ch3 cross section calculate: '
    real(musica_dk)    :: Tadj
    character(len=:), allocatable :: msg

    if( size(this%cross_section(1)%array,dim=2) == 4 ) then
      Tadj = min( 298._musica_dk,max( 235._musica_dk,environment%temperature ) )
      associate( coefficient => this%cross_section(1)%array )
      cross_section(:) = coefficient(:,1) &
                         *(rONE + Tadj*(coefficient(:,2) + Tadj*(coefficient(:,3) + Tadj*coefficient(:,4))))
      end associate
    else
      write(msg,*) Iam//'array must have 4 parameters'
      call die_msg( 500000001, msg )
    endif

    write(*,*) Iam,'exiting'

  end function run

end module micm_ch3coch3_ch3co_ch3_cross_section_type
