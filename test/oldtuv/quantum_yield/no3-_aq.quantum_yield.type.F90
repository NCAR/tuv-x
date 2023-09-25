! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: GPL-2.0-or-later
!
!> \file
!> This no3-aq+hv->no2(aq)+o- quantum yield module

!> The no3-aq+hv->no2(aq)+o- quantum yield type and related functions
module micm_no3m_aq_quantum_yield_type

  use micm_base_quantum_yield_type,    only : base_quantum_yield_t
  use musica_constants,                only : musica_dk, musica_ik

  implicit none

  private
  public :: no3m_aq_quantum_yield_t

  real(musica_dk), parameter :: rZERO = 0.0_musica_dk
  real(musica_dk), parameter :: rONE  = 1.0_musica_dk

  !> Calculator for no3m(aq)+hv->no2(aq)+o- quantum yield
  type, extends(base_quantum_yield_t) :: no3m_aq_quantum_yield_t
  contains
    !> Calculate the quantum yield
    procedure :: calculate => run
  end type no3m_aq_quantum_yield_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculate the photorate quantum yield for a given set of environmental conditions
  function run( this, environment ) result( quantum_yield )

    use micm_environment,                only : environment_t

    !> base quantum yield
    class(no3m_aq_quantum_yield_t), intent(in) :: this
    !> Environmental conditions
    class(environment_t), intent(in)                  :: environment
    !> Calculated quantum yield
    real(kind=musica_dk)                              :: quantum_yield(size(this%mdl_lambda_center))

    character(len=*), parameter :: Iam = 'no3-_(aq)+hv->products calculate: '

    write(*,*) Iam,'entering'

    quantum_yield = exp( -2400._musica_dk/environment%temperature + 3.6_musica_dk )      ! Chu & Anastasio, 2003

    write(*,*) Iam,'exiting'

  end function run

end module micm_no3m_aq_quantum_yield_type
