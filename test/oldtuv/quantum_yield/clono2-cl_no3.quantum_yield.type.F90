! Copyright (C) 2020-2025 University Corporation for Atmospheric Research
! SPDX-License-Identifier: GPL-2.0-or-later
!
!> \file
!> This clono2+hv->cl+no3 quantum yield module

!> The clono2+hv->cl+no3 quantum yield type and related functions
module micm_clono2_cl_no3_quantum_yield_type

  use micm_base_quantum_yield_type,    only : base_quantum_yield_t
  use musica_constants,                only : musica_dk, musica_ik

  implicit none

  private
  public :: clono2_cl_no3_quantum_yield_t

  real(musica_dk), parameter :: rZERO = 0.0_musica_dk
  real(musica_dk), parameter :: rONE  = 1.0_musica_dk

  !> Calculator for clono2+hv->oh+h quantum yield
  type, extends(base_quantum_yield_t) :: clono2_cl_no3_quantum_yield_t
  contains
    !> Calculate the quantum yield
    procedure :: calculate => run
  end type clono2_cl_no3_quantum_yield_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculate the photorate quantum yield for a given set of environmental conditions
  function run( this, environment ) result( quantum_yield )

    use micm_environment,                only : environment_t

    !> base quantum yield
    class(clono2_cl_no3_quantum_yield_t), intent(in) :: this
    !> Environmental conditions
    class(environment_t), intent(in)                  :: environment
    !> Calculated quantum yield
    real(kind=musica_dk)                              :: quantum_yield(size(this%mdl_lambda_center))

    character(len=*), parameter :: Iam = 'clono2+hv->cl+no3 calculate: '

    integer(musica_ik) :: wNdx
    real(musica_ik)    :: lambda

    do wNdx = 1,size(this%mdl_lambda_center)
      lambda = this%mdl_lambda_center(wNdx)
      if( lambda < 308._musica_dk ) then
        quantum_yield(wNdx) = .6_musica_dk
      elseif( 308._musica_dk <= lambda .and. lambda <= 364._musica_dk ) then
        quantum_yield(wNdx) = 7.143e-3_musica_dk * lambda - 1.6_musica_dk
      else
        quantum_yield(wNdx) = rONE
      endif
    enddo

  end function run

end module micm_clono2_cl_no3_quantum_yield_type
