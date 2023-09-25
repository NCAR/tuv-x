! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: GPL-2.0-or-later
!
!> \file
!> This n2o5+hv->no2_no3 cross_section module

!> The n2o5+hv->no2+no3_cross_section type and related functions
module micm_n2o5_no2_no3_cross_section_type

  use micm_base_cross_section_type,    only : base_cross_section_t

  implicit none

  private
  public :: n2o5_no2_no3_cross_section_t

  !> Calculator for base_cross_section
  type, extends(base_cross_section_t) :: n2o5_no2_no3_cross_section_t
  contains
    !> Calculate the cross section
    procedure :: calculate => run
  end type n2o5_no2_no3_cross_section_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculate the photorate cross section for a given set of environmental conditions
  function run( this, environment ) result( cross_section )

    use micm_environment, only : environment_t
    use musica_constants, only : musica_dk, musica_ik

    !> base cross section
    class(n2o5_no2_no3_cross_section_t), intent(in) :: this
    !> Environmental conditions
    class(environment_t), intent(in)                :: environment
    !> Calculated cross section
    real(kind=musica_dk)                            :: cross_section(size(this%mdl_lambda_center))

    character(len=*), parameter :: Iam = 'n2o5_no2_no3 cross section calculate: '
    real(musica_dk), parameter  :: Tlower = 233._musica_dk, Tupper = 300._musica_dk
    real(musica_dk), parameter  :: Tsf    = 1000._musica_dk
    integer(musica_ik) :: wNdx
    real(musica_dk) :: AdjTemp, Tfac

    write(*,*) Iam,'entering'

    associate( Temp => environment%temperature )
      AdjTemp = max( Tlower,min(Temp,Tupper) )
      do wNdx = 1,size(this%mdl_lambda_center)
        Tfac = Tsf*this%cross_section(2)%array(wNdx,1)*(Tupper - AdjTemp)/(Tupper*AdjTemp)
        cross_section(wNdx) = this%cross_section(1)%array(wNdx,1) * 10._musica_dk**(Tfac)
      enddo
    end associate

    write(*,*) Iam,'exiting'

  end function run

end module micm_n2o5_no2_no3_cross_section_type
