! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: GPL-2.0-or-later
!
!> \file
!> This ch3cocho+hv->ch3co+hco quantum yield module

!> The ch3cocho+hv->ch3co+hco quantum yield type and related functions
module micm_ch3cocho_ch3co_hco_quantum_yield_type

  use micm_base_quantum_yield_type,    only : base_quantum_yield_t
  use musica_constants,                only : musica_dk, musica_ik

  implicit none

  private
  public :: ch3cocho_ch3co_hco_quantum_yield_t

  real(musica_dk), parameter :: rZERO = 0.0_musica_dk
  real(musica_dk), parameter :: rONE  = 1.0_musica_dk

  !> Calculator for ch3cocho+hv->ch3co+hco quantum yield
  type, extends(base_quantum_yield_t) :: ch3cocho_ch3co_hco_quantum_yield_t
  contains
    !> Calculate the quantum yield
    procedure :: calculate => run
  end type ch3cocho_ch3co_hco_quantum_yield_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculate the photorate quantum yield for a given set of environmental conditions
  function run( this, environment ) result( quantum_yield )

    use micm_environment,                only : environment_t

    !> base quantum yield
    class(ch3cocho_ch3co_hco_quantum_yield_t), intent(in) :: this
    !> Environmental conditions
    class(environment_t), intent(in)                  :: environment
    !> Calculated quantum yield
    real(kind=musica_dk)                              :: quantum_yield(size(this%mdl_lambda_center))

    character(len=*), parameter :: Iam = 'ch3cocho+hv->ch3co_hco quantum yield calculate: '

    real(musica_dk), parameter :: lambdaL = 380._musica_dk
    real(musica_dk), parameter :: lambdaU = 440._musica_dk

    integer(musica_dk) :: wNdx
    real(musica_dk)    :: phi0, kq, lambda, airfac, qy
               
! zero pressure yield:
! 1.0 for wc < 380 nm
! 0.0 for wc > 440 nm
! linear in between:

! Pressure correction: quenching coefficient, torr-1
! in air, Koch and Moortgat:
    airfac = environment%number_density_air * 760._musica_dk/2.456E19_musica_dk

    do wNdx = 1,size(this%mdl_lambda_center)
      lambda = this%mdl_lambda_center(wNdx)
      phi0 = rONE - (lambda - 380._musica_dk)/60._musica_dk
      phi0 = MAX( MIN(phi0,rONE),rZERO )
      kq = 1.36e8_musica_dk * EXP( -8793._musica_dk/lambda )
      if( phi0 > rZERO ) then
        if( lambda >= lambdaL .and. lambda <= lambdaU ) then
          qy = phi0 / (phi0 + kq * airfac)
        else
          qy = phi0
        endif
      else
        qy = rZERO
      endif
      quantum_yield(wNdx) = qy
    enddo

    write(*,*) Iam,'exiting'

  end function run

end module micm_ch3cocho_ch3co_hco_quantum_yield_type
