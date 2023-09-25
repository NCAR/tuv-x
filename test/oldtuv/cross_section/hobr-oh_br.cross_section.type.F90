! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: GPL-2.0-or-later
!
!> \file
!> This hobr-oh_br cross_section module

!> The hobr-oh_br_cross_section type and related functions
module micm_hobr_oh_br_cross_section_type

  use micm_base_cross_section_type,    only : base_cross_section_t
  use musica_constants,                only : musica_dk, musica_ik

  implicit none

  private
  public :: hobr_oh_br_cross_section_t

  !> Calculator for hobr-oh_br cross section
  type, extends(base_cross_section_t) :: hobr_oh_br_cross_section_t
  contains
    !> Calculate the cross section
    procedure :: calculate => run
  end type hobr_oh_br_cross_section_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculate the photorate cross section for a given set of environmental conditions
  function run( this, environment ) result( cross_section )

    use micm_environment,                only : environment_t
    use musica_constants,                only : musica_dk

    !> hobr_oh_br cross section
    class(hobr_oh_br_cross_section_t), intent(in) :: this
    !> Environmental conditions
    class(environment_t), intent(in)               :: environment
    !> Calculated cross section
    real(kind=musica_dk)                           :: cross_section(size(this%mdl_lambda_center))

    real(musica_dk), parameter :: rZERO = 0.0_musica_dk
    real(musica_dk), parameter :: a = -2.359E-3_musica_dk
    real(musica_dk), parameter :: b = 1.2478_musica_dk
    real(musica_dk), parameter :: c = -210.4_musica_dk

    character(len=*), parameter :: Iam = 'hobr_oh_br cross section calculate: '

    write(*,*) Iam,'entering'

    associate( wc => this%mdl_lambda_center )
    where( wc >= 250._musica_dk .and. wc <= 550._musica_dk )
      cross_section = &
               24.77_musica_dk * exp( -109.80_musica_dk*(LOG(284.01_musica_dk/wc))**2 ) &
             + 12.22_musica_dk * exp(  -93.63_musica_dk*(LOG(350.57_musica_dk/wc))**2 ) &
             + 2.283_musica_dk * exp(- 242.40_musica_dk*(LOG(457.38_musica_dk/wc))**2 )
      cross_section = cross_section * 1.e-20_musica_dk
    elsewhere
      cross_section = rZERO
    endwhere
    end associate

    write(*,*) Iam,'exiting'

  end function run

end module micm_hobr_oh_br_cross_section_type
