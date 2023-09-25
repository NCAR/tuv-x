! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: GPL-2.0-or-later
!
!> \file
!> This nitroxy_acetone cross_section module

!> The nitroxy_acetone_cross_section type and related functions
module micm_nitroxy_acetone_cross_section_type

  use micm_base_cross_section_type,    only : base_cross_section_t
  use musica_constants,                only : musica_dk, musica_ik

  implicit none

  private
  public :: nitroxy_acetone_cross_section_t

  !> Calculator for nitroxy_acetone cross section
  type, extends(base_cross_section_t) :: nitroxy_acetone_cross_section_t
  contains
    !> Calculate the cross section
    procedure :: calculate => run
  end type nitroxy_acetone_cross_section_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculate the photorate cross section for a given set of environmental conditions
  function run( this, environment ) result( cross_section )

    use micm_environment,                only : environment_t
    use musica_constants,                only : musica_dk

    !> nitroxy_acetone cross section
    class(nitroxy_acetone_cross_section_t), intent(in) :: this
    !> Environmental conditions
    class(environment_t), intent(in)               :: environment
    !> Calculated cross section
    real(kind=musica_dk)                           :: cross_section(size(this%mdl_lambda_center))

    real(musica_dk), parameter :: rZERO = 0.0_musica_dk
    real(musica_dk), parameter :: a = -1.365E-3_musica_dk
    real(musica_dk), parameter :: b = 0.7834_musica_dk
    real(musica_dk), parameter :: c = -156.8_musica_dk

    character(len=*), parameter :: Iam = 'nitroxy_acetone cross section calculate: '

    write(*,*) Iam,'entering'

    where( this%mdl_lambda_center >= 284._musica_dk .and. this%mdl_lambda_center <= 335._musica_dk )
      cross_section = exp( c + this%mdl_lambda_center*(b + a*this%mdl_lambda_center) )
    elsewhere
      cross_section = rZERO
    endwhere

    write(*,*) Iam,'exiting'

  end function run

end module micm_nitroxy_acetone_cross_section_type
