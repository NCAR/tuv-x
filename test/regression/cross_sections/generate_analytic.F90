! Copyright (C) 2026 University Corporation for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
! Standalone Fortran reference data generator for the five analytic
! (temperature- and height-independent) cross-section algorithms.
! No TUV-x or NetCDF dependencies -- formulas extracted directly from
! main:src/cross_sections/{rayliegh,hobr-oh_br,t_butyl_nitrate,
!   nitroxy_acetone,nitroxy_ethanol}.F90
!
! Output: one CSV per species in the same directory as this program.
! Columns: wavelength_m, cross_section_m2
!
! Wavelength grid: 26 bins, 200--450 nm in 10 nm steps, covering all
! active ranges plus out-of-range bins for zero-weight verification.
!
! Compile and run:
!   gfortran -O2 -o generate_analytic generate_analytic.F90
!   ./generate_analytic

program generate_analytic

  implicit none

  integer, parameter :: dp = kind(1.0d0)

  ! 26 wavelength mid-points [nm], 200-450 nm in 10 nm steps
  integer, parameter :: nwl = 26
  real(dp), parameter :: wl_start_nm = 200.0_dp
  real(dp), parameter :: wl_step_nm  = 10.0_dp

  real(dp) :: wl_nm(nwl), wl_m(nwl)
  real(dp) :: sigma_cm2(nwl), sigma_m2(nwl)
  integer  :: i

  ! Build wavelength grid
  do i = 1, nwl
    wl_nm(i) = wl_start_nm + (i - 1) * wl_step_nm
    wl_m(i)  = wl_nm(i) * 1.0e-9_dp
  end do

  call write_rayleigh(wl_nm, wl_m, nwl)
  call write_hobr(wl_nm, wl_m, nwl)
  call write_nitrate(wl_nm, wl_m, nwl, &
       -0.993e-3_dp, 0.5307_dp, -115.5_dp, &
       270.0_dp, 330.0_dp, 't_butyl_nitrate.csv')
  call write_nitrate(wl_nm, wl_m, nwl, &
       -1.365e-3_dp, 0.7834_dp, -156.8_dp, &
       284.0_dp, 335.0_dp, 'nitroxy_acetone.csv')
  call write_nitrate(wl_nm, wl_m, nwl, &
       -2.359e-3_dp, 1.2478_dp, -210.4_dp, &
       270.0_dp, 306.0_dp, 'nitroxy_ethanol.csv')

  write(*,*) 'Reference data written.'

contains

  ! Rayleigh scattering (rayliegh.F90)
  ! wrk = lambda_nm / 1000  (in micrometers)
  ! pwr = 3.6772 + 0.389*wrk + 0.09426/wrk  (wrk <= 0.55)
  !     = 4.04                                (wrk > 0.55)
  ! sigma_cm2 = 4.02e-28 / wrk^pwr
  subroutine write_rayleigh(wl_nm, wl_m, n)
    real(dp), intent(in) :: wl_nm(n), wl_m(n)
    integer,  intent(in) :: n
    real(dp) :: wrk, pwr, sig_cm2, sig_m2
    integer  :: i, u
    open(newunit=u, file='rayleigh.csv', status='replace')
    write(u, '(a)') 'wavelength_m,cross_section_m2'
    do i = 1, n
      wrk = wl_nm(i) * 1.0e-3_dp
      if (wrk <= 0.55_dp) then
        pwr = 3.6772_dp + 0.389_dp * wrk + 0.09426_dp / wrk
      else
        pwr = 4.04_dp
      end if
      sig_cm2 = 4.02e-28_dp / (wrk ** pwr)
      sig_m2  = sig_cm2 * 1.0e-4_dp
      write(u, '(es24.17, a, es24.17)') wl_m(i), ',', sig_m2
    end do
    close(u)
  end subroutine write_rayleigh

  ! HOBr (hobr-oh_br.F90)
  ! sigma_cm2 = [24.77*exp(-109.80*(ln(284.01/wl))^2)
  !            + 12.22*exp( -93.63*(ln(350.57/wl))^2)
  !            + 2.283*exp(-242.40*(ln(457.38/wl))^2)] * 1e-20
  ! Active for 250 <= wl_nm <= 550, zero elsewhere
  subroutine write_hobr(wl_nm, wl_m, n)
    real(dp), intent(in) :: wl_nm(n), wl_m(n)
    integer,  intent(in) :: n
    real(dp) :: sig_cm2, sig_m2
    integer  :: i, u
    open(newunit=u, file='hobr.csv', status='replace')
    write(u, '(a)') 'wavelength_m,cross_section_m2'
    do i = 1, n
      if (wl_nm(i) >= 250.0_dp .and. wl_nm(i) <= 550.0_dp) then
        sig_cm2 = ( 24.77_dp  * exp(-109.80_dp * (log(284.01_dp / wl_nm(i)))**2) &
                  + 12.22_dp  * exp( -93.63_dp * (log(350.57_dp / wl_nm(i)))**2) &
                  + 2.283_dp  * exp(-242.40_dp * (log(457.38_dp / wl_nm(i)))**2) ) &
                  * 1.0e-20_dp
      else
        sig_cm2 = 0.0_dp
      end if
      sig_m2 = sig_cm2 * 1.0e-4_dp
      write(u, '(es24.17, a, es24.17)') wl_m(i), ',', sig_m2
    end do
    close(u)
  end subroutine write_hobr

  ! Generic nitrate/nitroxy formula (t_butyl_nitrate.F90, nitroxy_*.F90)
  ! sigma_cm2 = exp(c + wl_nm * (b + a * wl_nm))
  ! Active for wl_min_nm <= wl_nm <= wl_max_nm, zero elsewhere
  subroutine write_nitrate(wl_nm, wl_m, n, a, b, c, wl_min, wl_max, fname)
    real(dp), intent(in)      :: wl_nm(n), wl_m(n)
    integer,  intent(in)      :: n
    real(dp), intent(in)      :: a, b, c, wl_min, wl_max
    character(*), intent(in)  :: fname
    real(dp) :: sig_cm2, sig_m2
    integer  :: i, u
    open(newunit=u, file=fname, status='replace')
    write(u, '(a)') 'wavelength_m,cross_section_m2'
    do i = 1, n
      if (wl_nm(i) >= wl_min .and. wl_nm(i) <= wl_max) then
        sig_cm2 = exp(c + wl_nm(i) * (b + a * wl_nm(i)))
      else
        sig_cm2 = 0.0_dp
      end if
      sig_m2 = sig_cm2 * 1.0e-4_dp
      write(u, '(es24.17, a, es24.17)') wl_m(i), ',', sig_m2
    end do
    close(u)
  end subroutine write_nitrate

end program generate_analytic
