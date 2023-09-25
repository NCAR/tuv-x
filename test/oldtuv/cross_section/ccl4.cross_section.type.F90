! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: GPL-2.0-or-later
!
!> \file
!> This ccl4->Products cross_section module

!> The ccl4->Products_cross_section type and related functions
module micm_ccl4_cross_section_type

  use micm_base_cross_section_type,    only : base_cross_section_t
  use musica_constants,                only : musica_dk, musica_ik

  implicit none

  private
  public :: ccl4_cross_section_t

  !> Calculator for ccl4 cross section
  type, extends(base_cross_section_t) :: ccl4_cross_section_t
  contains
    !> Calculate the cross section
    procedure :: calculate => run
  end type ccl4_cross_section_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculate the photorate cross section for a given set of environmental conditions
  function run( this, environment ) result( cross_section )

    use micm_environment,                only : environment_t
    use musica_constants,                only : musica_dk

    !> ccl4->products cross section
    class(ccl4_cross_section_t), intent(in) :: this
    !> Environmental conditions
    class(environment_t), intent(in)        :: environment
    !> Calculated cross section
    real(kind=musica_dk)                    :: cross_section(size(this%mdl_lambda_center))

    real(musica_dk), parameter :: b0 = 1.0739_musica_dk
    real(musica_dk), parameter :: b1 = -1.6275e-2_musica_dk
    real(musica_dk), parameter :: b2 = 8.8141e-5_musica_dk
    real(musica_dk), parameter :: b3 = -1.9811e-7_musica_dk
    real(musica_dk), parameter :: b4 = 1.5022e-10_musica_dk

    character(len=*), parameter :: Iam = 'ccl4 cross section run: '
    integer(musica_ik) :: wNdx
    real(musica_dk)    :: Temp, lambda, Wpoly
    real(musica_dk)    :: w1, w2, w3, w4

    write(*,*) Iam,'entering'

    Temp = max( min( 300._musica_dk,environment%temperature ),210._musica_dk )
    Temp = Temp - 295._musica_dk 
    do wNdx = 1,size(this%mdl_lambda_center)
      w1 = this%mdl_lambda_center(wNdx)
      if( w1 > 194._musica_dk .and. w1 < 250._musica_dk ) then
        w2 = w1**2
        w3 = w1**3
        w4 = w1**4
        Wpoly = b0 + b1*w1 + b2*w2 + b3*w3 + b4*w4
        cross_section(wNdx) = this%cross_section(1)%array(wNdx,1) * 10._musica_dk**(Wpoly*Temp)
      else
        cross_section(wNdx) = this%cross_section(1)%array(wNdx,1)
      endif
    enddo

    write(*,*) Iam,'exiting'

  end function run

end module micm_ccl4_cross_section_type
