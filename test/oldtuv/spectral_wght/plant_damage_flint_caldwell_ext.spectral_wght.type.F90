! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: GPL-2.0-or-later
!
!> \file
!> This plant_damage_flint_caldwell_ext_spectral_wght module

!> The spectral_wght type and related functions
module micm_plant_damage_flint_caldwell_ext_spectral_wght_type

  use micm_base_spectral_wght_type,    only : base_spectral_wght_t
  use musica_constants,                only : musica_dk, musica_ik

  implicit none

  private
  public :: plant_damage_flint_caldwell_ext_spectral_wght_t

  !> Calculator for plant_damage_flint_caldwell_ext_spectral_wght
  type, extends(base_spectral_wght_t) :: plant_damage_flint_caldwell_ext_spectral_wght_t
  contains
    !> Calculate the spectral weight
    procedure :: calculate => run
  end type plant_damage_flint_caldwell_ext_spectral_wght_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculate the spectral weight for a given set of environmental conditions
  function run( this, environment ) result( spectral_wght )

    use micm_environment,                only : environment_t

    !> plant_damage_flint_caldwell_ext type
    class(plant_damage_flint_caldwell_ext_spectral_wght_t), intent(in) :: this
    !> Environmental conditions
    class(environment_t), intent(in)        :: environment
    !> Calculated spectral weight
    real(kind=musica_dk)                    :: spectral_wght(size(this%mdl_lambda_center))

    real(musica_dk), parameter  :: a0 = 4.688272_musica_dk
    real(musica_dk), parameter  :: a1 = .1703411_musica_dk
    real(musica_dk), parameter  :: w1 = 307.867_musica_dk
    real(musica_dk), parameter  :: w2 = 390._musica_dk
    character(len=*), parameter :: Iam = 'plant_damage_flint_caldwell_ext calculate: '

    write(*,*) Iam,'entering'

    spectral_wght = EXP( a0*EXP(-EXP(a1*(this%mdl_lambda_center - w1)/1.15_musica_dk)) &
                                     + ((w2 - this%mdl_lambda_center)/121.7557_musica_dk - 4.183832_musica_dk) )
    spectral_wght = spectral_wght * this%mdl_lambda_center / 300._musica_dk
    where( spectral_wght < 0.0_musica_dk .or. this%mdl_lambda_center > 390._musica_dk )
      spectral_wght = 0.0_musica_dk
    endwhere

    write(*,*) Iam,'exiting'

  end function run

end module micm_plant_damage_flint_caldwell_ext_spectral_wght_type
