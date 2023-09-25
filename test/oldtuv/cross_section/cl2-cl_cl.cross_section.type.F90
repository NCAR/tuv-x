! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: GPL-2.0-or-later
!
!> \file
!> This cl2_cl_cl_cross_section module

!> The cl2+hv->cl_cl cross_section type and related functions
module micm_cl2_cl_cl_cross_section_type

  use micm_base_cross_section_type,    only : base_cross_section_t

  implicit none

  private
  public :: cl2_cl_cl_cross_section_t

  !> Calculator for base_cross_section
  type, extends(base_cross_section_t) :: cl2_cl_cl_cross_section_t
  contains
    !> Calculate the cross section
    procedure :: calculate => run
  end type cl2_cl_cl_cross_section_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculate the photorate cross section for a given set of environmental conditions
  function run( this, environment ) result( cross_section )

    use micm_environment,  only : environment_t
    use musica_constants,  only : musica_dk, musica_ik

    !> base cross section
    class(cl2_cl_cl_cross_section_t), intent(in) :: this
    !> Environmental conditions
    class(environment_t), intent(in)             :: environment
    !> Calculated cross section
    real(kind=musica_dk)                         :: cross_section(size(this%mdl_lambda_center))

    real(musica_dk), parameter :: rONE = 1.0_musica_dk
    character(len=*), parameter :: Iam = 'cl2+hv->cl+cl cross section calculate: '
    integer(musica_ik) :: wNdx
    real(musica_dk) :: aa, bb, bbsq, alpha, ex1, ex2

    write(*,*) Iam,'entering'

    aa = 402.7_musica_dk/environment%temperature
    bb = exp( aa )
    bbsq  = bb * bb
    alpha = (bbsq - rONE)/(bbsq + rONE)

    associate( wc => this%mdl_lambda_center )
      do wNdx = 1,size(this%mdl_lambda_center)
        ex1 = 27.3_musica_dk  * exp(-99.0_musica_dk * alpha * (log(329.5_musica_dk/wc(wNdx)))**2)
        ex2 =  .932_musica_dk * exp(-91.5_musica_dk * alpha * (log(406.5_musica_dk/wc(wNdx)))**2)
        cross_section(wNdx) = 1.e-20_musica_dk * sqrt(alpha) * (ex1 + ex2)
      enddo
    end associate

    write(*,*) Iam,'exiting'

  end function run

end module micm_cl2_cl_cl_cross_section_type
