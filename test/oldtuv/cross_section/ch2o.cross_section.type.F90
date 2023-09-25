! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: GPL-2.0-or-later
!
!> \file
!> This ch2o cross_section module

!> The ch2o_cross_section type and related functions
module micm_ch2o_cross_section_type

  use micm_base_cross_section_type,    only : base_cross_section_t
  use musica_constants,                only : musica_dk, musica_ik

  implicit none

  private
  public :: ch2o_cross_section_t

  !> Calculator for acetone cross_section
  type, extends(base_cross_section_t) :: ch2o_cross_section_t
  contains
    !> Initialize the cross section
    procedure :: calculate => run
  end type ch2o_cross_section_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculate the photorate cross section for a given set of environmental conditions
  function run( this, environment ) result( cross_section )

    use micm_environment,  only : environment_t

    !> base cross section
    class(ch2o_cross_section_t), intent(in) :: this
    !> Environmental conditions
    class(environment_t), intent(in)        :: environment
    !> Calculated cross section
    real(kind=musica_dk)                    :: cross_section(size(this%mdl_lambda_center))

    character(len=*), parameter :: Iam = 'ch2o cross section calculate: '
    real(musica_dk)    :: Tadj

    Tadj = environment%temperature - 298._musica_dk
    cross_section = this%cross_section(1)%array(:,1) + this%cross_section(1)%array(:,2) * Tadj 

    write(*,*) Iam,'exiting'

  end function run

end module micm_ch2o_cross_section_type
