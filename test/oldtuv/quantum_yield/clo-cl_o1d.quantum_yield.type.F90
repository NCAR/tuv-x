! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: GPL-2.0-or-later
!
!> \file
!> This clo+hv->cl+o1d quantum yield module

!> The clo+hv->cl+o1d quantum yield type and related functions
module micm_clo_cl_o1d_quantum_yield_type

  use micm_base_quantum_yield_type,    only : base_quantum_yield_t
  use musica_constants,                only : musica_dk, musica_ik

  implicit none

  private
  public :: clo_cl_o1d_quantum_yield_t

  real(musica_dk), parameter :: rZERO = 0.0_musica_dk
  real(musica_dk), parameter :: rONE  = 1.0_musica_dk

  !> Calculator for clo+hv->cl+o1d quantum yield
  type, extends(base_quantum_yield_t) :: clo_cl_o1d_quantum_yield_t
  contains
    !> Calculate the quantum yield
    procedure :: calculate => run
  end type clo_cl_o1d_quantum_yield_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculate the photorate quantum yield for a given set of environmental conditions
  function run( this, environment ) result( quantum_yield )

    use micm_environment,                only : environment_t

    !> base quantum yield
    class(clo_cl_o1d_quantum_yield_t), intent(in) :: this
    !> Environmental conditions
    class(environment_t), intent(in)                  :: environment
    !> Calculated quantum yield
    real(kind=musica_dk)                              :: quantum_yield(size(this%mdl_lambda_center))

    character(len=*), parameter :: Iam = 'clo+hv->cl+o1d calculate: '

    write(*,*) Iam,'entering'

    where( this%mdl_lambda_center < 263.4_musica_dk )
      quantum_yield = rONE
    elsewhere
      quantum_yield = rZERO
    endwhere

    write(*,*) Iam,'exiting'

  end function run

end module micm_clo_cl_o1d_quantum_yield_type
