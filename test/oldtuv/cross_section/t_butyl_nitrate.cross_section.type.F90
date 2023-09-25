! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: GPL-2.0-or-later
!
!> \file
!> This t_butyl_nitrate cross_section module

!> The t_butyl_nitrate_cross_section type and related functions
module micm_t_butyl_nitrate_cross_section_type

  use micm_base_cross_section_type,    only : base_cross_section_t
  use musica_constants,                only : musica_dk, musica_ik

  implicit none

  private
  public :: t_butyl_nitrate_cross_section_t

  !> Calculator for t_butyl_nitrate cross section
  type, extends(base_cross_section_t) :: t_butyl_nitrate_cross_section_t
  contains
    !> Calculate the cross section
    procedure :: calculate => run
  end type t_butyl_nitrate_cross_section_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculate the photorate cross section for a given set of environmental conditions
  function run( this, environment ) result( cross_section )

    use micm_environment,                only : environment_t
    use musica_constants,                only : musica_dk

    !> t_butyl_nitrate cross section
    class(t_butyl_nitrate_cross_section_t), intent(in) :: this
    !> Environmental conditions
    class(environment_t), intent(in)               :: environment
    !> Calculated cross section
    real(kind=musica_dk)                           :: cross_section(size(this%mdl_lambda_center))

    real(musica_dk), parameter :: rZERO = 0.0_musica_dk
    real(musica_dk), parameter :: a = -0.993E-3_musica_dk
    real(musica_dk), parameter :: b = 0.5307_musica_dk
    real(musica_dk), parameter :: c = -115.5_musica_dk

    character(len=*), parameter :: Iam = 't_butyl_nitrate cross section calculate: '

    write(*,*) Iam,'entering'

    where( this%mdl_lambda_center >= 270._musica_dk .and. this%mdl_lambda_center <= 330._musica_dk )
      cross_section = exp( c + this%mdl_lambda_center*(b + a*this%mdl_lambda_center) )
    elsewhere
      cross_section = rZERO
    endwhere

    write(*,*) Iam,'exiting'

  end function run

end module micm_t_butyl_nitrate_cross_section_type
