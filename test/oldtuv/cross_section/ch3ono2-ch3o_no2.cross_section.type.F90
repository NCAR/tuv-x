! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: GPL-2.0-or-later
!
!> \file
!> This ch3ono2->ch3o+no2 cross_section module

!> The ch3ono2->ch3o+no2_cross_section type and related functions
module micm_ch3ono2_ch3o_no2_cross_section_type

  use micm_base_cross_section_type,    only : base_cross_section_t
  use musica_constants,                only : musica_dk, musica_ik

  implicit none

  private
  public :: ch3ono2_ch3o_no2_cross_section_t

  !> Calculator for ch3ono2-ch3o_no2 cross section
  type, extends(base_cross_section_t) :: ch3ono2_ch3o_no2_cross_section_t
  contains
    !> Calculate the cross section
    procedure :: calculate => run
  end type ch3ono2_ch3o_no2_cross_section_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculate the photorate cross section for a given set of environmental conditions
  function run( this, environment ) result( cross_section )

    use micm_environment,                only : environment_t
    use musica_constants,                only : musica_dk

    !> hno3->oh+no2 cross section
    class(ch3ono2_ch3o_no2_cross_section_t), intent(in) :: this
    !> Environmental conditions
    class(environment_t), intent(in)               :: environment
    !> Calculated cross section
    real(kind=musica_dk)                           :: cross_section(size(this%mdl_lambda_center))

    character(len=*), parameter :: Iam = 'ch3ono2->ch3o+no2 cross section run: '
    real(musica_dk) :: Temp

    write(*,*) Iam,'entering'

    Temp = environment%temperature - 298._musica_dk 
    cross_section = this%cross_section(1)%array(:,1)*exp( this%cross_section(1)%array(:,2)*Temp )

    write(*,*) Iam,'exiting'

  end function run

end module micm_ch3ono2_ch3o_no2_cross_section_type
