! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: GPL-2.0-or-later
!
!> \file
!> This gaussian_340_nm_10_nm_FWHM_spectral_wght module

!> The spectral_wght type and related functions
module micm_gaussian_340_nm_10_nm_FWHM_spectral_wght_type

  use micm_base_spectral_wght_type,    only : base_spectral_wght_t
  use musica_constants,                only : musica_dk, musica_ik

  implicit none

  private
  public :: gaussian_340_nm_10_nm_FWHM_spectral_wght_t

  !> Calculator for gaussian_340_nm_10_nm_FWHM_spectral_wght
  type, extends(base_spectral_wght_t) :: gaussian_340_nm_10_nm_FWHM_spectral_wght_t
  contains
    !> Calculate the spectral weight
    procedure :: calculate => run
  end type gaussian_340_nm_10_nm_FWHM_spectral_wght_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculate the spectral weight for a given set of environmental conditions
  function run( this, environment ) result( spectral_wght )

    use micm_environment,                only : environment_t

    !> gaussian_340_nm_10_nm_FWHM type
    class(gaussian_340_nm_10_nm_FWHM_spectral_wght_t), intent(in) :: this
    !> Environmental conditions
    class(environment_t), intent(in)        :: environment
    !> Calculated cross section
    real(kind=musica_dk)                    :: spectral_wght(size(this%mdl_lambda_center))

    real(kind=musica_dk)        :: accum
    character(len=*), parameter :: Iam = 'gaussian_340_nm_10_nm_FWHM calculate: '

    write(*,*) Iam,'entering'

    spectral_wght = exp( -(log(2._musica_dk)*.04_musica_dk*(this%mdl_lambda_center(:) - 340._musica_dk)**2) )
    accum = sum( spectral_wght )
    spectral_wght = spectral_wght/accum

    write(*,*) Iam,'exiting'

  end function run

end module micm_gaussian_340_nm_10_nm_FWHM_spectral_wght_type
