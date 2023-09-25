! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: GPL-2.0-or-later
!
!> \file
!> This phytoplankton_boucher_spectral_wght module

!> The spectral_wght type and related functions
module micm_phytoplankton_boucher_spectral_wght_type

  use micm_base_spectral_wght_type,    only : base_spectral_wght_t
  use musica_constants,                only : musica_dk, musica_ik

  implicit none

  private
  public :: phytoplankton_boucher_spectral_wght_t

  !> Calculator for phytoplankton_boucher_spectral_wght
  type, extends(base_spectral_wght_t) :: phytoplankton_boucher_spectral_wght_t
  contains
    !> Calculate the spectral weight
    procedure :: calculate => run
  end type phytoplankton_boucher_spectral_wght_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculate the spectral weight for a given set of environmental conditions
  function run( this, environment ) result( spectral_wght )

    use micm_environment,                only : environment_t

    !> phytoplankton_boucher type
    class(phytoplankton_boucher_spectral_wght_t), intent(in) :: this
    !> Environmental conditions
    class(environment_t), intent(in)        :: environment
    !> Calculated spectral weight
    real(kind=musica_dk)                    :: spectral_wght(size(this%mdl_lambda_center))

    real(musica_dk), parameter  :: em = -3.17e-6_musica_dk
    real(musica_dk), parameter  :: a  = 112.5_musica_dk
    real(musica_dk), parameter  :: b  = -.6223_musica_dk
    real(musica_dk), parameter  :: c  = 7.67e-4_musica_dk
    character(len=*), parameter :: Iam = 'phytoplankton_boucher calculate: '

    write(*,*) Iam,'entering'

    where( this%mdl_lambda_center > 290._musica_dk .and. this%mdl_lambda_center < 400._musica_dk )
      spectral_wght = em + exp( a + this%mdl_lambda_center*(b + this%mdl_lambda_center*c) )
    elsewhere
      spectral_wght = 0.0_musica_dk
    endwhere

    write(*,*) Iam,'exiting'

  end function run

end module micm_phytoplankton_boucher_spectral_wght_type
