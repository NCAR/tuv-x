! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: GPL-2.0-or-later
!
!> \file
!> This ho2+hv->oh+o quantum yield module

!> The ho2+hv->oh+h quantum yield type and related functions
module micm_ho2_oh_o_quantum_yield_type

  use micm_base_quantum_yield_type,    only : base_quantum_yield_t
  use musica_constants,                only : musica_dk, musica_ik

  implicit none

  private
  public :: ho2_oh_o_quantum_yield_t

  real(musica_dk), parameter :: rZERO = 0.0_musica_dk
  real(musica_dk), parameter :: rONE  = 1.0_musica_dk

  !> Calculator for ho2+hv->oh+h quantum yield
  type, extends(base_quantum_yield_t) :: ho2_oh_o_quantum_yield_t
  contains
    !> Calculate the quantum yield
    procedure :: calculate => run
  end type ho2_oh_o_quantum_yield_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculate the photorate quantum yield for a given set of environmental conditions
  function run( this, environment ) result( quantum_yield )

    use micm_environment,                only : environment_t

    !> base quantum yield
    class(ho2_oh_o_quantum_yield_t), intent(in) :: this
    !> Environmental conditions
    class(environment_t), intent(in)                  :: environment
    !> Calculated quantum yield
    real(kind=musica_dk)                              :: quantum_yield(size(this%mdl_lambda_center))

    real(musica_dk), parameter :: lambda0 = 193._musica_dk
    character(len=*), parameter :: Iam = 'ho2+hv->oh+o calculate: '

    write(*,*) Iam,'entering'

    where( this%mdl_lambda_center >= 248._musica_dk )
      quantum_yield = rONE
    elsewhere
      quantum_yield = (rONE + 14._musica_dk*(this%mdl_lambda_center - lambda0)/55._musica_dk)/15._musica_dk
      quantum_yield = max( rZERO,quantum_yield )
    endwhere

    write(*,*) Iam,'exiting'

  end function run

end module micm_ho2_oh_o_quantum_yield_type
