! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: GPL-2.0-or-later
!
!> \file
!> This n2o+hv->n2+o1d cross_section module

!> The n2o+hv->n2+o1d_cross_section type and related functions
module micm_n2o_n2_o1d_cross_section_type

  use micm_base_cross_section_type,    only : base_cross_section_t

  implicit none

  private
  public :: n2o_n2_o1d_cross_section_t

  !> Calculator for base_cross_section
  type, extends(base_cross_section_t) :: n2o_n2_o1d_cross_section_t
  contains
    !> Calculate the cross section
    procedure :: calculate => run
  end type n2o_n2_o1d_cross_section_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculate the photorate cross section for a given set of environmental conditions
  function run( this, environment ) result( cross_section )

    use micm_environment, only : environment_t
    use musica_constants, only : musica_dk, musica_ik

    !> base cross section
    class(n2o_n2_o1d_cross_section_t), intent(in) :: this
    !> Environmental conditions
    class(environment_t), intent(in)                :: environment
    !> Calculated cross section
    real(kind=musica_dk)                            :: cross_section(size(this%mdl_lambda_center))

    character(len=*), parameter :: Iam = 'n2o_n2+o1d cross section calculate: '

    real(musica_dk), parameter :: rZERO = 0._musica_dk
    real(musica_dk), parameter :: A0 = 68.21023_musica_dk
    real(musica_dk), parameter :: A1 = -4.071805_musica_dk
    real(musica_dk), parameter :: A2 = 4.301146E-02_musica_dk
    real(musica_dk), parameter :: A3 = -1.777846E-04_musica_dk
    real(musica_dk), parameter :: A4 = 2.520672E-07_musica_dk

    real(musica_dk), parameter :: B0 = 123.4014_musica_dk
    real(musica_dk), parameter :: B1 = -2.116255_musica_dk
    real(musica_dk), parameter :: B2 = 1.111572E-02_musica_dk
    real(musica_dk), parameter :: B3 = -1.881058E-05_musica_dk

    integer(musica_ik) :: wNdx
    real(musica_dk) :: lambda, Tadj, A, B

!*** quantum yield of N(4s) and NO(2Pi) is less than 1% (Greenblatt and
!*** Ravishankara), so quantum yield of O(1D) is assumed to be unity

    Tadj = MAX(194._musica_dk,MIN(environment%temperature,320._musica_dk))
    do wNdx = 1,size(this%mdl_lambda_center)
      lambda = this%mdl_lambda_center(wNdx)   
      if( lambda >= 173._musica_dk .AND. lambda <= 240._musica_dk) then
        A = (((A4*lambda+A3)*lambda+A2)*lambda+A1)*lambda+A0
        B = (((B3*lambda+B2)*lambda+B1)*lambda+B0)
        B = (Tadj - 300._musica_dk)*EXP(B)
        cross_section(wNdx) = EXP(A+B)
      else
        cross_section(wNdx) = rZERO
      endif
    enddo

    write(*,*) Iam,'exiting'

  end function run

end module micm_n2o_n2_o1d_cross_section_type
