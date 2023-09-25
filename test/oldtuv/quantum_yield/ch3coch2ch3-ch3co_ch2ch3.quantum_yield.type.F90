! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: GPL-2.0-or-later
!
!> \file
!> This ch3coch2ch3+hv->ch3co+ch2ch3 quantum yield module

!> The ch3coch2ch3+hv->ch3co+ch2ch3 quantum yield type and related functions
module micm_ch3coch2ch3_quantum_yield_type

  use micm_base_quantum_yield_type,    only : base_quantum_yield_t
  use musica_constants,                only : musica_dk, musica_ik

  implicit none

  private
  public :: ch3coch2ch3_quantum_yield_t

  real(musica_dk), parameter :: rZERO = 0.0_musica_dk
  real(musica_dk), parameter :: rONE  = 1.0_musica_dk

  !> Calculator for ch3coch2ch3+hv->ch3co+ch2ch3 quantum yield
  type, extends(base_quantum_yield_t) :: ch3coch2ch3_quantum_yield_t
  contains
    !> Calculate the quantum yield
    procedure :: calculate => run
  end type ch3coch2ch3_quantum_yield_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculate the photorate quantum yield for a given set of environmental conditions
  function run( this, environment ) result( quantum_yield )

    use micm_environment,                only : environment_t

    !> base quantum yield
    class(ch3coch2ch3_quantum_yield_t), intent(in) :: this
    !> Environmental conditions
    class(environment_t), intent(in)                  :: environment
    !> Calculated quantum yield
    real(kind=musica_dk)                              :: quantum_yield(size(this%mdl_lambda_center))

    character(len=*), parameter :: Iam = 'ch3coch2ch3+hv->ch3co+ch2ch3 calculate: '
    real(musica_dk) :: ptorr

    write(*,*) Iam,'entering'

    ptorr = 760._musica_dk*environment%number_density_air/2.69e19_musica_dk
    quantum_yield = rONE/(0.96_musica_dk + 2.22E-3_musica_dk*ptorr)
    quantum_yield = min(quantum_yield, rONE)

    write(*,*) Iam,'exiting'

  end function run

end module micm_ch3coch2ch3_quantum_yield_type
