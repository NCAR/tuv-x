! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: GPL-2.0-or-later
!
!> \file
!> This clono2+hv cross_section module

!> The clono2+hv->n2+o1d_cross_section type and related functions
module micm_clono2_cross_section_type

  use micm_base_cross_section_type,    only : base_cross_section_t

  implicit none

  private
  public :: clono2_cross_section_t

  !> Calculator for clono2_cross_section
  type, extends(base_cross_section_t) :: clono2_cross_section_t
  contains
    !> Calculate the cross section
    procedure :: calculate => run
  end type clono2_cross_section_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculate the photorate cross section for a given set of environmental conditions
  function run( this, environment ) result( cross_section )

    use micm_environment, only : environment_t
    use musica_constants, only : musica_dk, musica_ik

    !> base cross section
    class(clono2_cross_section_t), intent(in) :: this
    !> Environmental conditions
    class(environment_t), intent(in)                :: environment
    !> Calculated cross section
    real(kind=musica_dk)                            :: cross_section(size(this%mdl_lambda_center))

    character(len=*), parameter :: Iam = 'clono2 cross section calculate: '

    real(musica_dk), parameter :: rONE = 1.0_musica_dk

    integer(musica_ik) :: wNdx
    real(musica_dk) :: Tadj

    write(*,*) Iam,'entering'

    Tadj = environment%temperature - 296._musica_dk

    associate( polyCoeff => this%cross_section(1)%array )
      cross_section = polyCoeff(:,1)*(rONE + Tadj*(polyCoeff(:,2) + Tadj*polyCoeff(:,3)))
    end associate

    write(*,*) Iam,'exiting'

  end function run

end module micm_clono2_cross_section_type
