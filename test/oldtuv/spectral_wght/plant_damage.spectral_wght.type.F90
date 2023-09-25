! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: GPL-2.0-or-later
!
!> \file
!> This plant_damage_spectral_wght module

!> The spectral_wght type and related functions
module micm_plant_damage_spectral_wght_type

  use micm_base_spectral_wght_type,    only : base_spectral_wght_t
  use musica_constants,                only : musica_dk, musica_ik

  implicit none

  private
  public :: plant_damage_spectral_wght_t

  !> Calculator for plant_damage_spectral_wght
  type, extends(base_spectral_wght_t) :: plant_damage_spectral_wght_t
  contains
    !> Calculate the spectral weight
    procedure :: calculate => run
  end type plant_damage_spectral_wght_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculate the spectral weight for a given set of environmental conditions
  function run( this, environment ) result( spectral_wght )

    use micm_environment,                only : environment_t

    !> plant_damage type
    class(plant_damage_spectral_wght_t), intent(in) :: this
    !> Environmental conditions
    class(environment_t), intent(in)        :: environment
    !> Calculated spectral weight
    real(kind=musica_dk)                    :: spectral_wght(size(this%mdl_lambda_center))

    real(musica_dk), parameter  :: a0 = 570.25_musica_dk
    real(musica_dk), parameter  :: a1 = -4.70144_musica_dk
    real(musica_dk), parameter  :: a2 = .01274_musica_dk
    real(musica_dk), parameter  :: a3 = -1.13118e-5_musica_dk
    character(len=*), parameter :: Iam = 'plant_damage calculate: '

    write(*,*) Iam,'entering'

    spectral_wght = a0 + this%mdl_lambda_center*(a1 + this%mdl_lambda_center*(a2 + this%mdl_lambda_center*a3))
    where( spectral_wght < 0.0_musica_dk .or. this%mdl_lambda_center > 313._musica_dk )
      spectral_wght = 0.0_musica_dk
    endwhere

    write(*,*) Iam,'exiting'

  end function run

end module micm_plant_damage_spectral_wght_type
