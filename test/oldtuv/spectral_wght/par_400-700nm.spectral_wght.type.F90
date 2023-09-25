! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: GPL-2.0-or-later
!
!> \file
!> This par_400_700nm_spectral_wght module

!> The spectral_wght type and related functions
module micm_par_400_700nm_spectral_wght_type

  use micm_base_spectral_wght_type,    only : base_spectral_wght_t
  use musica_constants,                only : musica_dk, musica_ik

  implicit none

  private
  public :: par_400_700nm_spectral_wght_t

  !> Calculator for par_400-700nm_spectral_wght
  type, extends(base_spectral_wght_t) :: par_400_700nm_spectral_wght_t
    !> The spectral wght array
  contains
    !> Calculate the spectral weight
    procedure :: calculate => run
  end type par_400_700nm_spectral_wght_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculate the spectral weight for a given set of environmental conditions
  function run( this, environment ) result( spectral_wght )

    use micm_environment,                only : environment_t

    !> par_400_700nm type
    class(par_400_700nm_spectral_wght_t), intent(in) :: this
    !> Environmental conditions
    class(environment_t), intent(in)        :: environment
    !> Calculated spectral weight
    real(kind=musica_dk)                    :: spectral_wght(size(this%mdl_lambda_center))

    character(len=*), parameter :: Iam = 'par_400_700nm calculate: '

    write(*,*) Iam,'entering'

    where( 400._musica_dk < this%mdl_lambda_center(:) .and. this%mdl_lambda_center(:) < 700._musica_dk )
      spectral_wght = 8.36e-3_musica_dk*this%mdl_lambda_center
    elsewhere
      spectral_wght = 0.0_musica_dk
    endwhere

    write(*,*) Iam,'exiting'

  end function run

end module micm_par_400_700nm_spectral_wght_type
