! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: GPL-2.0-or-later
!
!> \file
!> This clono2+hv->cl+no3 quantum yield module

!> The clono2+hv->cl+no3 quantum yield type and related functions
module micm_clono2_clo_no2_quantum_yield_type

  use micm_base_quantum_yield_type,    only : base_quantum_yield_t
  use musica_constants,                only : musica_dk, musica_ik

  implicit none

  private
  public :: clono2_clo_no2_quantum_yield_t

  real(musica_dk), parameter :: rZERO = 0.0_musica_dk
  real(musica_dk), parameter :: rONE  = 1.0_musica_dk

  !> Calculator for clono2+hv->clo+no2 quantum yield
  type, extends(base_quantum_yield_t) :: clono2_clo_no2_quantum_yield_t
  contains
    !> Calculate the quantum yield
    procedure :: calculate => run
  end type clono2_clo_no2_quantum_yield_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculate the photorate quantum yield for a given set of environmental conditions
  function run( this, environment ) result( quantum_yield )

    use micm_environment,                only : environment_t

    !> base quantum yield
    class(clono2_clo_no2_quantum_yield_t), intent(in) :: this
    !> Environmental conditions
    class(environment_t), intent(in)                  :: environment
    !> Calculated quantum yield
    real(kind=musica_dk)                              :: quantum_yield(size(this%mdl_lambda_center))

    character(len=*), parameter :: Iam = 'clono2+hv->clo+no2 calculate: '

    integer(musica_ik) :: wNdx
    real(musica_ik)    :: lambda, qyield

    write(*,*) Iam,'entering'

    do wNdx = 1,size(this%mdl_lambda_center)
      lambda = this%mdl_lambda_center(wNdx)
      if( lambda < 308._musica_dk ) then
        qyield = .6_musica_dk
      elseif( 308._musica_dk <= lambda .and. lambda <= 364._musica_dk ) then
        qyield = 7.143e-3_musica_dk * lambda - 1.6_musica_dk
      else
        qyield = rONE
      endif
      quantum_yield(wNdx) = rONE - qyield
    enddo

    write(*,*) Iam,'exiting'

  end function run

end module micm_clono2_clo_no2_quantum_yield_type
