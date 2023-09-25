! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: GPL-2.0-or-later
!
!> \file
!> This visplus_spectral_wght module

!> The spectral_wght type and related functions
module micm_visplus_spectral_wght_type

  use micm_base_spectral_wght_type,    only : base_spectral_wght_t
  use musica_constants,                only : musica_dk, musica_ik

  implicit none

  private
  public :: visplus_spectral_wght_t

  !> Calculator for uv-b_280-320_nm_spectral_wght
  type, extends(base_spectral_wght_t) :: visplus_spectral_wght_t
    !> The spectral wght array
  contains
    !> Calculate the cross section
    procedure :: calculate => run
  end type visplus_spectral_wght_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculate the photorate cross section for a given set of environmental conditions
  function run( this, environment ) result( spectral_wght )

    use micm_environment,                only : environment_t

    !> ub_b_315_400_nm type
    class(visplus_spectral_wght_t), intent(in) :: this
    !> Environmental conditions
    class(environment_t), intent(in)        :: environment
    !> Calculated cross section
    real(kind=musica_dk)                    :: spectral_wght(size(this%mdl_lambda_center))

    character(len=*), parameter :: Iam = 'visplus calculate: '

    write(*,*) Iam,'entering'

    where( 400._musica_dk < this%mdl_lambda_center(:) .and. this%mdl_lambda_center(:) < 700._musica_dk )
      spectral_wght = 1.0_musica_dk
    elsewhere
      spectral_wght = 0.0_musica_dk
    endwhere

    write(*,*) Iam,'exiting'

  end function run

end module micm_visplus_spectral_wght_type
