! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: GPL-2.0-or-later
!
!> \file
!> This cfc-11->Products cross_section module

!> The cfc-11->Products type and related functions
module micm_cfc11_cross_section_type

  use micm_base_cross_section_type,    only : base_cross_section_t
  use musica_constants,                only : musica_dk, musica_ik

  implicit none

  private
  public :: cfc11_cross_section_t

  !> Calculator for cfc-11->Products cross section
  type, extends(base_cross_section_t) :: cfc11_cross_section_t
  contains
    !> Calculate the cross section
    procedure :: calculate => run
  end type cfc11_cross_section_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculate the photorate cross section for a given set of environmental conditions
  function run( this, environment ) result( cross_section )

    use micm_environment,                only : environment_t
    use musica_constants,                only : musica_dk

    !> hno3->oh+no2 cross section
    class(cfc11_cross_section_t), intent(in) :: this
    !> Environmental conditions
    class(environment_t), intent(in)               :: environment
    !> Calculated cross section
    real(kind=musica_dk)                           :: cross_section(size(this%mdl_lambda_center))

    character(len=*), parameter :: Iam = 'cfc-11->Products cross section run: '
    real(musica_dk) :: Temp

    write(*,*) Iam,'entering'

    Temp = 1.e-4_musica_dk*(environment%temperature - 298._musica_dk)
    cross_section = this%cross_section(1)%array(:,1)*exp( (this%mdl_lambda_center(:) - 184.9_musica_dk)*Temp )

    write(*,*) Iam,'exiting'

  end function run

end module micm_cfc11_cross_section_type
