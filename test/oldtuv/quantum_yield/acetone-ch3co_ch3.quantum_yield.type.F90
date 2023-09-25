! Copyright (C) 2020 National Center for Atmospheric Research
! SPDX-License-Identifier: GPL-2.0-or-later
!
!> \file
!> This ch3coch3+hv->ch3co_ch3 quantum_yield module

!> The ch3coch3+hv->ch3co+ch3_quantum_yield type and related functions
module micm_ch3coch3_ch3co_ch3_quantum_yield_type

  use micm_base_quantum_yield_type,    only : base_quantum_yield_t
  use musica_constants,                only : musica_dk, musica_ik

  implicit none

  private
  public :: ch3coch3_ch3co_ch3_quantum_yield_t

  !> Calculator for acetone quantum_yield
  type, extends(base_quantum_yield_t) :: ch3coch3_ch3co_ch3_quantum_yield_t
  contains
    !> Initialize the quantum_yield
    procedure :: calculate => run
  end type ch3coch3_ch3co_ch3_quantum_yield_t

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !> Calculate the photorate quantum_yield for a given set of environmental conditions
  function run( this, environment ) result( quantum_yield )
!     qyacet - q.y. for acetone, based on Blitz et al. (2004)
! Compute acetone quantum yields according to the parameterization of:
! Blitz, M. A., D. E. Heard, M. J. Pilling, S. R. Arnold, and M. P. Chipperfield 
!       (2004), Pressure and temperature-dependent quantum yields for the 
!       photodissociation of acetone between 279 and 327.5 nm, Geophys. 
!       Res. Lett., 31, L06111, doi:10.1029/2003GL018793.

    use micm_environment,  only : environment_t

    !> base quantum_yield
    class(ch3coch3_ch3co_ch3_quantum_yield_t), intent(in) :: this
    !> Environmental conditions
    class(environment_t), intent(in)             :: environment
    !> Calculated quantum_yield
    real(kind=musica_dk)                         :: quantum_yield(size(this%mdl_lambda_center))

    real(musica_dk), parameter :: rZERO = 0.0_musica_dk
    real(musica_dk), parameter :: rONE  = 1.0_musica_dk
    character(len=*), parameter :: Iam = 'ch3coch3+hv->ch3co+ch3 quantum_yield calculate: '
    integer(musica_ik) :: wNdx
! w = wavelength, nm
! T = temperature, K
! M = air number density, molec. cm-3
    real(musica_dk)    :: w, wadj, Tadj, M
    real(musica_dk)    :: a0, a1, a2, a3, a4
    real(musica_dk)    :: b0, b1, b2, b3, b4
    real(musica_dk)    :: c3
    real(musica_dk)    :: cA0, cA1, cA2, cA3, cA4
    real(musica_dk)    :: dumexp
    real(musica_dk)    :: fco, fac

    write(*,*) Iam,'entering'

    Tadj = environment%temperature/295._musica_dk
    M    = environment%number_density_air
lambda_loop: &
    do wNdx = 1,size(this%mdl_lambda_center)
      w = this%mdl_lambda_center(wNdx)
      if(w < 279._musica_dk) then
         fac = 0.95_musica_dk
      elseif(w > 327._musica_dk) then
         fac = rZERO
      else
!   CO (carbon monoxide) quantum yields:
        a0 = 0.350_musica_dk * Tadj**(-1.28_musica_dk)
        b0 = 0.068_musica_dk * Tadj**(-2.65_musica_dk)
! SM: prevent exponent overflow in rare cases:
        dumexp = b0*(w - 248._musica_dk)
        if (dumexp > 80._musica_dk) then
           cA0 = 5.e34_musica_dk
        else
           cA0 = exp(dumexp) * a0 / (rONE - a0)
        endif

        fco = rONE / (rONE + cA0)
!   CH3CO (acetyl radical) quantum yields:
        wadj = 1.e7_musica_dk/w
        if(w >= 279._musica_dk .and. w < 302._musica_dk) then
          a1 = 1.600E-19_musica_dk * Tadj**(-2.38_musica_dk)
          b1 = 0.55E-3_musica_dk   * Tadj**(-3.19_musica_dk)
          cA1 = a1 * EXP(-b1*(wadj - 33113._musica_dk))
          fac = (rONE - fco) / (rONE + cA1 * M)
        elseif(w >= 302._musica_dk .AND. w <= 327._musica_dk) then
          a2 = 1.62E-17_musica_dk * Tadj**(-10.03_musica_dk)
          b2 = 1.79E-3_musica_dk  * Tadj**(-1.364_musica_dk)
          cA2 = a2 * EXP(-b2*(wadj - 30488._musica_dk))

          a3 = 26.29_musica_dk   * Tadj**(-6.59_musica_dk)
          b3 = 5.72E-7_musica_dk * Tadj**(-2.93_musica_dk)
          c3 = 30006._musica_dk  * Tadj**(-0.064_musica_dk)
          ca3 = a3 * EXP(-b3*(wadj - c3)**2)

          a4 = 1.67E-15_musica_dk * Tadj**(-7.25_musica_dk)
          b4 = 2.08E-3_musica_dk  * Tadj**(-1.16_musica_dk)
          cA4 = a4 * EXP(-b4*(wadj - 30488._musica_dk))

          fac = (rONE - fco) * (rONE + cA3 + cA4 * M) &
                / ((rONE + cA3 + cA2 * M)*(rONE + cA4 * M))

        endif
      endif
      quantum_yield(wNdx) = fac
    enddo lambda_loop

    write(*,*) Iam,'exiting'

  end function run

end module micm_ch3coch3_ch3co_ch3_quantum_yield_type
