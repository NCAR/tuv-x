! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: GPL-2.0-or-later
!
!> \file
!> This chcl3+hv->products cross section module

!> The chcl3+hv->products cross section type and related functions
module micm_chcl3_cross_section_type

  use micm_base_cross_section_type,    only : base_cross_section_t
  use musica_constants,                only : musica_dk, musica_ik

  implicit none

  private
  public :: chcl3_cross_section_t

  real(musica_dk), parameter :: rZERO = 0.0_musica_dk
  real(musica_dk), parameter :: rONE  = 1.0_musica_dk

  !> Calculator for chcl3+hv->oh+h cross section
  type, extends(base_cross_section_t) :: chcl3_cross_section_t
  contains
    !> Calculate the cross section
    procedure :: calculate => run
  end type chcl3_cross_section_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculate the photorate cross section for a given set of environmental conditions
  function run( this, environment ) result( cross_section )

    use micm_environment,                only : environment_t

    !> base cross section
    class(chcl3_cross_section_t), intent(in) :: this
    !> Environmental conditions
    class(environment_t), intent(in)                  :: environment
    !> Calculated cross section
    real(kind=musica_dk)                              :: cross_section(size(this%mdl_lambda_center))

    real(musica_dk), parameter :: b0 = 3.7973_musica_dk
    real(musica_dk), parameter :: b1 = -7.0913e-2_musica_dk
    real(musica_dk), parameter :: b2 = 4.9397e-4_musica_dk
    real(musica_dk), parameter :: b3 = -1.5226e-6_musica_dk
    real(musica_dk), parameter :: b4 = 1.7555e-9_musica_dk

    character(len=*), parameter :: Iam = 'chcl3+hv->products calculate: '
    integer(musica_ik) :: wndx
    real(musica_dk) :: w1, w2, w3, w4, tcoeff, Tadj

    write(*,*) Iam,'entering'

    associate( wc => this%mdl_lambda_center, Temp => environment%temperature )
      Tadj = min(max(Temp,210._musica_dk),300._musica_dk) - 295._musica_dk
      do wNdx = 1,size(wc)
        cross_section(wNdx) = this%cross_section(1)%array(wNdx,1)
        if( wc(wNdx) > 190._musica_dk .and. wc(wNdx) < 240._musica_dk ) then
          w1 = wc(wNdx)
!         w2 = w1**2
!         w3 = w1**3
!         w4 = w1**4
!         tcoeff = b0 + b1*w1 + b2*w2 + b3*w3 + b4*w4
          tcoeff = b0 + w1*(b1 + w1*(b2 + w1*(b3 + w1*b4)))
          cross_section(wNdx) = cross_section(wNdx) * 10._musica_dk**(tcoeff*Tadj)
        endif
      enddo
    end associate

    write(*,*) Iam,'exiting'

  end function run

end module micm_chcl3_cross_section_type
