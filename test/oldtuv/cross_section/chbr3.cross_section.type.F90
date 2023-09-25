! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: GPL-2.0-or-later
!
!> \file
!> This chbr3_cross_section module

!> The chbr3+hv->products cross_section type and related functions
module micm_chbr3_cross_section_type

  use micm_base_cross_section_type,    only : base_cross_section_t

  implicit none

  private
  public :: chbr3_cross_section_t

  !> Calculator for base_cross_section
  type, extends(base_cross_section_t) :: chbr3_cross_section_t
  contains
    !> Calculate the cross section
    procedure :: calculate => run
  end type chbr3_cross_section_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculate the photorate cross section for a given set of environmental conditions
  function run( this, environment ) result( cross_section )

    use micm_environment,  only : environment_t
    use musica_constants,  only : musica_dk, musica_ik

    !> base cross section
    class(chbr3_cross_section_t), intent(in) :: this
    !> Environmental conditions
    class(environment_t), intent(in)             :: environment
    !> Calculated cross section
    real(kind=musica_dk)                         :: cross_section(size(this%mdl_lambda_center))

    real(musica_dk), parameter :: C0 = .06183_musica_dk
    real(musica_dk), parameter :: C1 = .000241_musica_dk
    real(musica_dk), parameter :: C2 = 2.376_musica_dk
    real(musica_dk), parameter :: C3 = .14757_musica_dk
    real(musica_dk), parameter :: T0 = 273._musica_dk

    character(len=*), parameter :: Iam = 'chbr3+hv->Products cross section calculate: '
    integer(musica_ik) :: wNdx
    real(musica_dk)    :: wc, Temp, lambda

    write(*,*) Iam,'entering'

    do wNdx = 1,size(this%mdl_lambda_center)
      Temp = environment%temperature
      lambda = this%mdl_lambda_center(wNdx)
      if( lambda > 290._musica_dk .and. lambda < 340._musica_dk .and. &
          Temp > 210._musica_dk .and. Temp < 300._musica_dk ) then
        cross_section(wNdx) = &
          EXP( (C0 - C1*lambda)*(T0 - Temp) - (C2 + C3*lambda) )
      else
        cross_section(wNdx) = this%cross_section(1)%array(wNdx,1)
      endif
    enddo

    write(*,*) Iam,'exiting'

  end function run

end module micm_chbr3_cross_section_type
