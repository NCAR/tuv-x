! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: GPL-2.0-or-later
!
!> \file
!> This exponential_decay_spectral_wght module

!> The spectral_wght type and related functions
module micm_exponential_decay_spectral_wght_type

  use micm_base_spectral_wght_type,    only : base_spectral_wght_t
  use musica_constants,                only : musica_dk, musica_ik

  implicit none

  private
  public :: exponential_decay_spectral_wght_t

  !> Calculator for exponential_decay_spectral_wght
  type, extends(base_spectral_wght_t) :: exponential_decay_spectral_wght_t
    !> The spectral wght array
  contains
    !> Calculate the spectral weight
    procedure :: calculate => run
  end type exponential_decay_spectral_wght_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculate the spectral weight for a given set of environmental conditions
  function run( this, environment ) result( spectral_wght )

    use micm_environment,                only : environment_t

    !> exponential_decay type
    class(exponential_decay_spectral_wght_t), intent(in) :: this
    !> Environmental conditions
    class(environment_t), intent(in)        :: environment
    !> Calculated spectral weight
    real(kind=musica_dk)                    :: spectral_wght(size(this%mdl_lambda_center))

    character(len=*), parameter :: Iam = 'exponential_decay calculate: '

    write(*,*) Iam,'entering'

    spectral_wght = 10._musica_dk**((300._musica_dk - this%mdl_lambda_center)/14._musica_dk)

    write(*,*) Iam,'exiting'

  end function run

end module micm_exponential_decay_spectral_wght_type
