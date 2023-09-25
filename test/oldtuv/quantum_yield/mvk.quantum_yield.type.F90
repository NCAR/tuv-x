! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: GPL-2.0-or-later
!
!> \file
!> This mvk+hv->products quantum yield module

!> The mvk+hv->prodcuts quantum yield type and related functions
module micm_mvk_quantum_yield_type

  use micm_base_quantum_yield_type,    only : base_quantum_yield_t
  use musica_constants,                only : musica_dk, musica_ik

  implicit none

  private
  public :: mvk_quantum_yield_t

  real(musica_dk), parameter :: rZERO = 0.0_musica_dk
  real(musica_dk), parameter :: rONE  = 1.0_musica_dk

  !> Calculator for mvk+hv->oh+h quantum yield
  type, extends(base_quantum_yield_t) :: mvk_quantum_yield_t
  contains
    !> Calculate the quantum yield
    procedure :: calculate => run
  end type mvk_quantum_yield_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculate the photorate quantum yield for a given set of environmental conditions
  function run( this, environment ) result( quantum_yield )

    use micm_environment,                only : environment_t

    !> base quantum yield
    class(mvk_quantum_yield_t), intent(in) :: this
    !> Environmental conditions
    class(environment_t), intent(in)                  :: environment
    !> Calculated quantum yield
    real(kind=musica_dk)                              :: quantum_yield(size(this%mdl_lambda_center))

    character(len=*), parameter :: Iam = 'mvk+hv->products calculate: '
    real(musica_dk) :: divisor

    write(*,*) Iam,'entering'

    divisor = 5.5_musica_dk + 9.2e-19_musica_dk*environment%number_density_air
    quantum_yield = exp( -0.055_musica_dk*(this%mdl_lambda_center - 308._musica_dk)) / divisor
    quantum_yield = min( quantum_yield,rONE )

    write(*,*) Iam,'exiting'

  end function run

end module micm_mvk_quantum_yield_type
