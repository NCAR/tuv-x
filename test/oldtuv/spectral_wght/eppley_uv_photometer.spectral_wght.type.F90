! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: GPL-2.0-or-later
!
!> \file
!> This eppley_uv_photometer_spectral_wght module

!> The spectral_wght type and related functions
module micm_eppley_uv_photometer_spectral_wght_type

  use micm_base_spectral_wght_type,    only : base_spectral_wght_t
  use musica_constants,                only : musica_dk, musica_ik

  implicit none

  private
  public :: eppley_uv_photometer_spectral_wght_t

  !> Calculator for eppley_uv_photometer_spectral_wght
  type, extends(base_spectral_wght_t) :: eppley_uv_photometer_spectral_wght_t
  contains
    !> Calculate the spectral wght
    procedure :: calculate => run
  end type eppley_uv_photometer_spectral_wght_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculate the spectral wght for a given set of environmental conditions
  function run( this, environment ) result( spectral_wght )

    use micm_environment,                only : environment_t

    !> eppley_uv_photometer spectral wght
    class(eppley_uv_photometer_spectral_wght_t), intent(in) :: this
    !> Environmental conditions
    class(environment_t), intent(in)        :: environment
    !> Calculated spectral wght
    real(kind=musica_dk)                    :: spectral_wght(size(this%mdl_lambda_center))

    integer                     :: nLambda
    real(kind=musica_dk)        :: accum
    character(len=*), parameter :: Iam = 'eppley_uv_photometer spectral wght calculate: '

    write(*,*) Iam,'entering'

    nLambda = size(this%mdl_lambda_edge)
    spectral_wght = this%spectral_wght(1)%array(:,1)
    accum = sum( spectral_wght*(this%mdl_lambda_edge(2:nLambda) - this%mdl_lambda_edge(1:nLambda-1)) )
    spectral_wght = 90._musica_dk*spectral_wght/accum

    write(*,*) Iam,'exiting'

  end function run

end module micm_eppley_uv_photometer_spectral_wght_type
