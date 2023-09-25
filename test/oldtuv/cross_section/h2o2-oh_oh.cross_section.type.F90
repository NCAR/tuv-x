! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: GPL-2.0-or-later
!
!> \file
!> This h2o2_oh_oh_cross_section module

!> The h2o2+hv->oh_oh cross_section type and related functions
module micm_h2o2_oh_oh_cross_section_type

  use micm_base_cross_section_type,    only : base_cross_section_t

  implicit none

  private
  public :: h2o2_oh_oh_cross_section_t

  !> Calculator for base_cross_section
  type, extends(base_cross_section_t) :: h2o2_oh_oh_cross_section_t
  contains
    !> Calculate the cross section
    procedure :: calculate => run
  end type h2o2_oh_oh_cross_section_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculate the photorate cross section for a given set of environmental conditions
  function run( this, environment ) result( cross_section )

    use micm_environment,  only : environment_t
    use musica_constants,  only : musica_dk, musica_ik

    !> base cross section
    class(h2o2_oh_oh_cross_section_t), intent(in) :: this
    !> Environmental conditions
    class(environment_t), intent(in)             :: environment
    !> Calculated cross section
    real(kind=musica_dk)                         :: cross_section(size(this%mdl_lambda_center))

    real(musica_dk), parameter :: rONE = 1.0_musica_dk
    real(musica_dk), parameter ::  A0 = 6.4761E+04_musica_dk
    real(musica_dk), parameter ::  A1 = -9.2170972E+02_musica_dk
    real(musica_dk), parameter ::  A2 = 4.535649_musica_dk
    real(musica_dk), parameter ::  A3 = -4.4589016E-03_musica_dk
    real(musica_dk), parameter ::  A4 = -4.035101E-05_musica_dk
    real(musica_dk), parameter ::  A5 = 1.6878206E-07_musica_dk
    real(musica_dk), parameter ::  A6 = -2.652014E-10_musica_dk
    real(musica_dk), parameter ::  A7 = 1.5534675E-13_musica_dk

    real(musica_dk), parameter ::  B0 = 6.8123E+03_musica_dk
    real(musica_dk), parameter ::  B1 = -5.1351E+01_musica_dk
    real(musica_dk), parameter ::  B2 = 1.1522E-01_musica_dk
    real(musica_dk), parameter ::  B3 = -3.0493E-05_musica_dk
    real(musica_dk), parameter ::  B4 = -1.0924E-07_musica_dk

    character(len=*), parameter :: Iam = 'h2o2+hv->oh+oh cross section calculate: '
    integer(musica_ik) :: wNdx
    real(musica_dk) :: lambda, sumA, sumB, t, chi, xs

    write(*,*) Iam,'entering'

    associate( wl => this%mdl_lambda_edge, wc => this%mdl_lambda_center )
      do wNdx = 1,size(this%mdl_lambda_center)
! Parameterization (JPL94)
! Range 260-350 nm; 200-400 K
        if( wl(wNdx) >= 260._musica_dk .and. wl(wNdx) < 350._musica_dk ) then
           lambda = wc(wNdx)
           sumA = ((((((A7*lambda + A6)*lambda + A5)*lambda + A4)*lambda +A3)*lambda + A2)*lambda + A1)*lambda + A0
           sumB = (((B4*lambda + B3)*lambda + B2)*lambda + B1)*lambda + B0
           t = MIN(MAX(environment%temperature,200._musica_dk),400._musica_dk)            
           chi = rONE/(rONE + EXP(-1265._musica_dk/t))
           cross_section(wNdx) = (chi * sumA + (rONE - chi)*sumB)*1.E-21_musica_dk
         else
           cross_section(wNdx) = this%cross_section(1)%array(wNdx,1)
         endif
      enddo
    end associate

    write(*,*) Iam,'exiting'

  end function run

end module micm_h2o2_oh_oh_cross_section_type
