! Copyright (C) 2023-2026 University Corporation for Atmospheric Research
! SPDX-License-Identifier: Apache-2.0
!
! Standalone Fortran reference-data generator for delta-Eddington regression tests.
!
! Extracts the core algorithm from:
!   src/radiative_transfer/solvers/delta_eddington.F90
!   src/spherical_geometry.F90
!   src/linear_algebras/linpack.F90
! without any framework dependencies (no musica, no warehouses).
!
! Test atmosphere: 3-layer, tau=0.5/layer, ssa=0.9/layer, g=0.85/layer,
!                 albedo=0.1, altitude edges [0,1,2,3] km.
!
! Usage:
!   gfortran -O2 -o generate_reference generate_reference.F90
!   ./generate_reference
!
! Output: delta_eddington_sza_X.XXXX.csv (one per SZA in sza_list)

program generate_reference

  implicit none

  integer,  parameter :: dp       = kind(1.0d0)
  real(dp), parameter :: pi       = 3.1415926535898_dp
  real(dp), parameter :: d2r      = pi / 180.0_dp
  real(dp), parameter :: re       = 6.371e3_dp   ! Earth radius [km]
  real(dp), parameter :: largest  = 1.0e36_dp
  real(dp), parameter :: rZERO    = 0.0_dp
  real(dp), parameter :: rONE     = 1.0_dp
  real(dp), parameter :: rTWO     = 2.0_dp
  real(dp), parameter :: NINETY   = 90.0_dp

  integer,  parameter :: n_layers = 3
  integer,  parameter :: n_levels = n_layers + 1

  ! Test atmosphere (uniform across all layers; indexing top-to-bottom, layer 1 = topmost)
  real(dp) :: tauu(n_layers), omu_in(n_layers), gu_in(n_layers)
  real(dp) :: rsfc

  ! Altitude edges bottom-to-top [km]: 0, 1, 2, 3
  real(dp) :: ze(n_levels)      ! bottom-to-top
  real(dp) :: zd(0:n_layers)    ! top-to-bottom (zd(0)=TOA, zd(n_layers)=ground)

  ! Spherical-geometry output
  integer  :: nid(0:n_layers)
  real(dp) :: dsdh(0:n_layers, n_layers)

  ! Radiation field (indexed 1=TOA at entry, reversed to 1=ground at output)
  real(dp) :: edr(n_levels), eup(n_levels), edn(n_levels)
  real(dp) :: fdr(n_levels), fup(n_levels), fdn(n_levels)

  ! SZAs to process [radians]
  real(dp) :: sza_list(3)

  integer          :: isza, lev, unit
  real(dp)         :: sza_rad, zen_deg, mu
  character(len=256) :: fname

  ! --- Set test atmosphere ---
  tauu   = 0.5_dp
  omu_in = 0.9_dp
  gu_in  = 0.85_dp
  rsfc   = 0.1_dp

  ! Altitude edges [km], bottom-to-top: 0, 1, 2, 3
  ze(1) = 0.0_dp; ze(2) = 1.0_dp; ze(3) = 2.0_dp; ze(4) = 3.0_dp
  ! Reverse to top-to-bottom: zd(0)=3, zd(1)=2, zd(2)=1, zd(3)=0
  do lev = 0, n_layers
    zd(lev) = ze(n_levels - lev)
  end do

  sza_list = [0.0_dp, 0.5_dp, 1.0_dp]

  do isza = 1, 3
    sza_rad = sza_list(isza)
    zen_deg = sza_rad / d2r
    mu      = cos(sza_rad)

    call set_params(zen_deg, n_layers, zd, nid, dsdh)
    call solve_delta_eddington(n_layers, mu, rsfc, tauu, omu_in, gu_in, &
                               nid, dsdh, edr, eup, edn, fdr, fup, fdn)

    write(fname, '("delta_eddington_sza_",f6.4,".csv")') sza_rad
    open(newunit=unit, file=trim(fname), status='replace', action='write')
    write(unit,'(a)') '"# atmosphere: 3-layer, tau=0.5/layer, ssa=0.9/layer,' &
                    //' g=0.85/layer, albedo=0.1"'
    write(unit,'("# sza_rad=",f3.1)') sza_rad
    write(unit,'(a)') '"# level 0 = ground, level n_layers = TOA"'
    write(unit,'(a)') 'level,direct_irradiance,upwelling_irradiance,' &
                    //'downwelling_irradiance,direct_actinic_flux,' &
                    //'upwelling_actinic_flux,downwelling_actinic_flux'
    ! After reversal in solve_delta_eddington: index 1=ground (CSV level 0),
    !                                          index n_levels=TOA (CSV level n_layers)
    do lev = 0, n_layers
      write(unit,'(i0,",",es22.15e2,",",es22.15e2,",",es22.15e2,",",'// &
                         'es22.15e2,",",es22.15e2,",",es22.15e2)') &
            lev, edr(lev+1), eup(lev+1), edn(lev+1), &
                 fdr(lev+1), fup(lev+1), fdn(lev+1)
    end do
    close(unit)
    write(*,'("Written: ",a)') trim(fname)
  end do

contains

  !---------------------------------------------------------------------------
  ! Adapted from spherical_geometry.F90 :: set_parameters
  ! Computes nid(0:n_lay) and dsdh(0:n_lay, n_lay) for the direct-beam path.
  !---------------------------------------------------------------------------
  subroutine set_params(zen_deg, n_lay, zd_in, nid_out, dsdh_out)

    real(dp), intent(in)  :: zen_deg
    integer,  intent(in)  :: n_lay
    real(dp), intent(in)  :: zd_in(0:n_lay)   ! altitudes top-to-bottom [km]
    integer,  intent(out) :: nid_out(0:n_lay)
    real(dp), intent(out) :: dsdh_out(0:n_lay, n_lay)

    real(dp) :: zenrad, sinrad, rpsinz
    real(dp) :: rj, rjp1, dhj, ga, gb, dsj, sm
    integer  :: i, j, id

    zenrad = zen_deg * d2r
    sinrad = sin(zenrad)
    nid_out   = 0
    dsdh_out  = rZERO

    do i = 0, n_lay
      rpsinz = (re + zd_in(i)) * sinrad

      if (zen_deg > NINETY .and. rpsinz < re) then
        id = -1
      else
        id = i
        if (zen_deg > NINETY) then
          id = -1
          do j = 1, n_lay
            if (rpsinz < (zd_in(j-1) + re) .and. &
                rpsinz >= (zd_in(j)   + re)) id = j
          end do
        end if

        do j = 1, id
          sm = rONE
          if (j == id .and. id == i .and. zen_deg > NINETY) sm = -rONE
          rj   = re + zd_in(j-1)
          rjp1 = re + zd_in(j)
          dhj  = zd_in(j-1) - zd_in(j)
          ga = rj*rj - rpsinz*rpsinz
          gb = rjp1*rjp1 - rpsinz*rpsinz
          ga = max(rZERO, ga)
          gb = max(rZERO, gb)
          if (id > i .and. j == id) then
            dsj = sqrt(ga)
          else
            dsj = sqrt(ga) - sm * sqrt(gb)
          end if
          dsdh_out(i, j) = dsj / dhj
        end do
      end if

      nid_out(i) = id
    end do

  end subroutine set_params

  !---------------------------------------------------------------------------
  ! Adapted from solver.F90 :: slant_optical_depth
  ! Returns 1e36 when nid_i < 0 (below tangent height).
  !---------------------------------------------------------------------------
  pure real(dp) function slant_od(i_lay, nid_i, dsdh_row, taun)

    integer,  intent(in) :: i_lay
    integer,  intent(in) :: nid_i
    real(dp), intent(in) :: dsdh_row(:)
    real(dp), intent(in) :: taun(:)

    integer :: j

    if (nid_i < 0) then
      slant_od = largest
      return
    end if
    slant_od = rZERO
    do j = 1, min(nid_i, i_lay)
      slant_od = slant_od + taun(j) * dsdh_row(j)
    end do
    do j = min(nid_i, i_lay) + 1, nid_i
      slant_od = slant_od + 2.0_dp * taun(j) * dsdh_row(j)
    end do

  end function slant_od

  !---------------------------------------------------------------------------
  ! Adapted from linear_algebras/linpack.F90 :: tridiag
  ! Solves B*u = r where B has lower diagonal a, main diagonal b, upper c.
  !---------------------------------------------------------------------------
  subroutine tridiag_solve(a, b, c, r, u, n)

    integer,  intent(in)  :: n
    real(dp), intent(in)  :: a(n), b(n), c(n), r(n)
    real(dp), intent(out) :: u(n)

    integer  :: i
    real(dp) :: denom
    real(dp) :: cp(n)

    cp(1) = c(1) / b(1)
    u(1)  = r(1) / b(1)
    do i = 2, n
      denom = rONE / (b(i) - a(i) * cp(i-1))
      cp(i) = c(i) * denom
      u(i)  = (r(i) - a(i) * u(i-1)) * denom
    end do
    do i = n - 1, 1, -1
      u(i) = u(i) - cp(i) * u(i+1)
    end do

  end subroutine tridiag_solve

  !---------------------------------------------------------------------------
  ! Adapted from delta_eddington.F90 :: update_radiation_field (inner kernel)
  ! Single wavelength, single column.
  ! Inputs:  n_lay, mu, rsfc, tauu, omu_arr, gu_arr, nid, dsdh
  ! Outputs: edr/eup/edn (spectral irradiance), fdr/fup/fdn (actinic flux)
  !          Arrays indexed 1=ground, n_lay+1=TOA (bottom-to-top)
  !---------------------------------------------------------------------------
  subroutine solve_delta_eddington(n_lay, mu, rsfc_in, tauu_in, omu_arr, &
                                   gu_arr, nid_in, dsdh_in,              &
                                   edr, eup, edn, fdr, fup, fdn)

    integer,  intent(in)  :: n_lay
    real(dp), intent(in)  :: mu, rsfc_in
    real(dp), intent(in)  :: tauu_in(n_lay), omu_arr(n_lay), gu_arr(n_lay)
    integer,  intent(in)  :: nid_in(0:n_lay)
    real(dp), intent(in)  :: dsdh_in(0:n_lay, n_lay)
    real(dp), intent(out) :: edr(n_lay+1), eup(n_lay+1), edn(n_lay+1)
    real(dp), intent(out) :: fdr(n_lay+1), fup(n_lay+1), fdn(n_lay+1)

    real(dp), parameter :: precis = 1.0e-7_dp
    real(dp), parameter :: eps    = 1.0e-3_dp

    real(dp) :: tausla(0:n_lay), tauc(0:n_lay), mu2(0:n_lay)
    real(dp) :: gi(n_lay), omi(n_lay), taun(n_lay)
    real(dp) :: lam(n_lay), bgam(n_lay), mu1(n_lay)
    real(dp) :: e1(n_lay), e2(n_lay), e3(n_lay), e4(n_lay)
    real(dp) :: cup(n_lay), cdn(n_lay), cuptn(n_lay), cdntn(n_lay)

    integer  :: mrows
    real(dp), allocatable :: la(:), lb(:), lc(:), rhs(:), y(:)

    real(dp) :: pifs, fdn0, surfem, ssfc
    real(dp) :: f, g, om, tempg
    real(dp) :: gam1, gam2, gam3, gam4
    real(dp) :: expon, expon0, expon1, divisr, temp, up, dn
    integer  :: i, j, row, lev

    pifs   = rONE
    fdn0   = rZERO
    surfem = rZERO

    tauc   = rZERO
    tausla = rZERO
    mu2    = rONE / sqrt(largest)

    ! Delta-scaling (Eddington approximation)
    do i = 1, n_lay
      f       = gu_arr(i) * gu_arr(i)
      gi(i)   = (gu_arr(i) - f)   / (rONE - f)
      omi(i)  = (rONE - f) * omu_arr(i) / (rONE - omu_arr(i) * f)
      taun(i) = (rONE - omu_arr(i) * f) * tauu_in(i)
    end do

    ! TOA slant depth for below-horizon case
    if (mu < rZERO) then
      tausla(0) = slant_od(0, nid_in(0), dsdh_in(0,:), taun)
    end if

    do i = 1, n_lay
      g  = gi(i)
      om = omi(i)
      tauc(i) = tauc(i-1) + taun(i)

      tempg = min(abs(g), rONE - precis)
      g = sign(tempg, g)
      om = min(om, rONE - precis)

      tausla(i) = slant_od(i, nid_in(i), dsdh_in(i,:), taun)

      if (nid_in(i) >= 0) then
        if (tausla(i) == tausla(i-1)) then
          mu2(i) = sqrt(largest)
        else
          mu2(i) = (tauc(i) - tauc(i-1)) / (tausla(i) - tausla(i-1))
          mu2(i) = sign(max(abs(mu2(i)), rONE / sqrt(largest)), mu2(i))
        end if
      end if

      ! Eddington coefficients (Toon et al. 1989, Table 1)
      gam1 =   (7.0_dp - om * (4.0_dp + 3.0_dp * g)) / 4.0_dp
      gam2 = - (rONE   - om * (4.0_dp - 3.0_dp * g)) / 4.0_dp
      gam3 =   (rTWO   - 3.0_dp * g * mu) / 4.0_dp
      gam4 =   rONE - gam3
      mu1(i) = 0.5_dp

      lam(i) = sqrt(gam1*gam1 - gam2*gam2)

      if (gam2 /= rZERO) then
        bgam(i) = (gam1 - lam(i)) / gam2
      else
        bgam(i) = rZERO
      end if

      expon = exp(-lam(i) * taun(i))
      e1(i) = rONE + bgam(i) * expon
      e2(i) = rONE - bgam(i) * expon
      e3(i) = bgam(i) + expon
      e4(i) = bgam(i) - expon

      expon0 = exp(-tausla(i-1))
      expon1 = exp(-tausla(i))

      divisr = lam(i)*lam(i) - rONE/(mu2(i)*mu2(i))
      temp   = max(eps, abs(divisr))
      divisr = sign(temp, divisr)

      up = om * pifs * ((gam1 - rONE/mu2(i)) * gam3 + gam4 * gam2) / divisr
      dn = om * pifs * ((gam1 + rONE/mu2(i)) * gam4 + gam2 * gam3) / divisr

      cup(i)   = up * expon0
      cdn(i)   = dn * expon0
      cuptn(i) = up * expon1
      cdntn(i) = dn * expon1
    end do

    ssfc  = rsfc_in * mu * exp(-tausla(n_lay)) * pifs + surfem
    mrows = 2 * n_lay

    allocate(la(mrows), lb(mrows), lc(mrows), rhs(mrows), y(mrows))

    ! First row
    la(1)  = rZERO
    lb(1)  = e1(1)
    lc(1)  = -e2(1)
    rhs(1) = fdn0 - cdn(1)

    ! Odd rows 3 .. mrows-1
    i = 0
    do row = 3, mrows-1, 2
      i = i + 1
      la(row)  = e2(i) * e3(i)   - e4(i) * e1(i)
      lb(row)  = e1(i) * e1(i+1) - e3(i) * e3(i+1)
      lc(row)  = e3(i) * e4(i+1) - e1(i) * e2(i+1)
      rhs(row) = e3(i) * (cup(i+1) - cuptn(i)) + e1(i) * (cdntn(i) - cdn(i+1))
    end do

    ! Even rows 2 .. mrows-2
    i = 0
    do row = 2, mrows-2, 2
      i = i + 1
      la(row)  = e2(i+1) * e1(i)   - e3(i)   * e4(i+1)
      lb(row)  = e2(i)   * e2(i+1) - e4(i)   * e4(i+1)
      lc(row)  = e1(i+1) * e4(i+1) - e2(i+1) * e3(i+1)
      rhs(row) = (cup(i+1) - cuptn(i)) * e2(i+1) &
               - (cdn(i+1) - cdntn(i)) * e4(i+1)
    end do

    ! Last row
    la(mrows)  = e1(n_lay) - rsfc_in * e3(n_lay)
    lb(mrows)  = e2(n_lay) - rsfc_in * e4(n_lay)
    lc(mrows)  = rZERO
    rhs(mrows) = ssfc - cuptn(n_lay) + rsfc_in * cdntn(n_lay)

    call tridiag_solve(la, lb, lc, rhs, y, mrows)

    ! Unfold (top-to-bottom: level 1=TOA, level n_lay+1=ground)
    fdr(1) = pifs * exp(-tausla(0))
    edr(1) = mu * fdr(1)
    edn(1) = fdn0
    eup(1) = y(1) * e3(1) - y(2) * e4(1) + cup(1)
    fdn(1) = edn(1) / mu1(1)
    fup(1) = eup(1) / mu1(1)

    j   = 1
    row = 1
    do lev = 2, n_lay + 1
      fdr(lev) = pifs * exp(-tausla(lev-1))
      edr(lev) = mu * fdr(lev)
      edn(lev) = y(row) * e3(j) + y(row+1) * e4(j) + cdntn(j)
      eup(lev) = y(row) * e1(j) + y(row+1) * e2(j) + cuptn(j)
      fdn(lev) = edn(lev) / mu1(j)
      fup(lev) = eup(lev) / mu1(j)
      row = row + 2
      j   = j + 1
    end do

    ! Reverse to bottom-to-top (index 1=ground, n_lay+1=TOA)
    edr = edr(n_lay+1:1:-1)
    eup = eup(n_lay+1:1:-1)
    edn = edn(n_lay+1:1:-1)
    fdr = fdr(n_lay+1:1:-1)
    fup = fup(n_lay+1:1:-1)
    fdn = fdn(n_lay+1:1:-1)

    deallocate(la, lb, lc, rhs, y)

  end subroutine solve_delta_eddington

end program generate_reference
