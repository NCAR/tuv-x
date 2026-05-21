#!/usr/bin/env python3
# Copyright (C) 2023-2026 University Corporation for Atmospheric Research
# SPDX-License-Identifier: Apache-2.0
#
# Generates reference CSV data for the delta-Eddington regression tests.
#
# This is an independent Python port of Toon et al. (1989) following the same
# algorithm as the C++ implementation.  Run this script whenever the test
# atmosphere or SZA list changes, commit the resulting CSV files, and the C++
# regression test will compare against them automatically.
#
# Reference:
#   Toon et al., J. Geophys. Res. 94(D13), 16287-16301, 1989.
#   doi:10.1029/JD094iD13p16287
#   Dahlback & Stamnes, Planet. Space Sci. 39(5), 671-683, 1991.
#   doi:10.1016/0032-0633(91)90061-E

import csv
import math
import os

# ── Test atmosphere ────────────────────────────────────────────────────────────

# 3-layer atmosphere with 1 wavelength band, 1 column
# All optical properties are top-to-bottom (index 0 = topmost layer)
ALT_EDGES_M   = [0.0, 1000.0, 2000.0, 3000.0]   # bottom-to-top [m]
OPTICAL_DEPTH = [0.5, 0.5, 0.5]                  # per layer, top-to-bottom
SINGLE_SCATTERING_ALBEDO = [0.9, 0.9, 0.9]
ASYMMETRY_PARAMETER      = [0.85, 0.85, 0.85]
SURFACE_ALBEDO = 0.1

# SZAs to generate reference data for [radians]
SZAS = [0.0, 0.5, 1.0]

# ── Spherical geometry ─────────────────────────────────────────────────────────

EARTH_RADIUS_M = 6.371e6


def _compute_geometry(sza_rad, alt_edges_m):
    """Return (nid, dsdh) matching SphericalGeometry::SetParameters."""
    n_layers = len(alt_edges_m) - 1
    n_levels = n_layers + 1
    surface_elevation = alt_edges_m[0]
    re = EARTH_RADIUS_M + surface_elevation
    half_pi = math.pi / 2.0
    above_horizon = sza_rad <= half_pi
    sin_sza = math.sin(sza_rad)

    # Reorder: altitude_from_toa[0] = TOA height, [n_layers] = 0
    altitude_from_toa = [alt_edges_m[n_layers - k] - surface_elevation for k in range(n_levels)]

    nid  = [0] * n_levels
    dsdh = [[0.0] * n_layers for _ in range(n_levels)]

    for i in range(n_levels):
        impact_parameter = (re + altitude_from_toa[i]) * sin_sza

        if not above_horizon and impact_parameter < re:
            layers_crossed = -1
        else:
            layers_crossed = i
            if not above_horizon:
                layers_crossed = -1
                for j in range(1, n_layers + 1):
                    if (impact_parameter < altitude_from_toa[j - 1] + re and
                            impact_parameter >= altitude_from_toa[j] + re):
                        layers_crossed = j

            for j in range(1, layers_crossed + 1):
                path_sign    = 1.0
                radius_upper = re + altitude_from_toa[j - 1]
                radius_lower = re + altitude_from_toa[j]
                layer_thickness = altitude_from_toa[j - 1] - altitude_from_toa[j]
                upper_term = max(0.0, radius_upper**2 - impact_parameter**2)
                lower_term = max(0.0, radius_lower**2 - impact_parameter**2)

                if j == layers_crossed == i and not above_horizon:
                    path_sign = -1.0

                slant_length = path_sign * (math.sqrt(upper_term) - math.sqrt(lower_term))
                if layer_thickness > 0.0:
                    dsdh[i][j - 1] = slant_length / layer_thickness

        nid[i] = layers_crossed

    return nid, dsdh


def _slant_optical_depth(level, n_layers_crossed, slant_path, optical_depth):
    """Return slant optical depth matching SlantOpticalDepth()."""
    if n_layers_crossed < 0:
        return math.inf
    result = 0.0
    direct_path_layers = min(n_layers_crossed, level)
    for j in range(direct_path_layers):
        result += optical_depth[j] * slant_path[j]
    for j in range(direct_path_layers, n_layers_crossed):
        result += 2.0 * optical_depth[j] * slant_path[j]
    return result


# ── Delta-Eddington solver ─────────────────────────────────────────────────────

def _thomas_solve(lower, main, upper, rhs):
    """Thomas algorithm matching tuvx::Solve(TridiagonalMatrix, Array1D)."""
    n = len(rhs)
    main  = list(main)
    rhs   = list(rhs)
    for i in range(1, n):
        factor = lower[i - 1] / main[i - 1]
        main[i] -= factor * upper[i - 1]
        rhs[i]  -= factor * rhs[i - 1]
    y = [0.0] * n
    y[-1] = rhs[-1] / main[-1]
    for i in range(n - 2, -1, -1):
        y[i] = (rhs[i] - upper[i] * y[i + 1]) / main[i]
    return y


def solve(sza_rad, optical_depth, single_scattering_albedo, asymmetry_parameter,
          surface_albedo, alt_edges_m):
    """
    Run the delta-Eddington solver for a single column and wavelength.

    Parameters — all top-to-bottom (index 0 = topmost layer) except alt_edges_m
    which is bottom-to-top (index 0 = ground).

    Returns a dict keyed by field name, each value a list of length n_levels
    ordered bottom-to-top (index 0 = ground, index n_layers = TOA).
    """
    LARGEST   = 1.0e36
    PRECISION = 1.0e-7
    EPSILON   = 1.0e-3

    n_layers = len(optical_depth)
    n_levels = n_layers + 1
    mu = math.cos(sza_rad)

    nid, dsdh = _compute_geometry(sza_rad, alt_edges_m)

    # Delta-scaling
    scaled_asymmetry = [0.0] * n_layers
    scaled_ssa       = [0.0] * n_layers
    scaled_tau       = [0.0] * n_layers
    for i in range(n_layers):
        g     = asymmetry_parameter[i]
        omega = single_scattering_albedo[i]
        tau   = optical_depth[i]
        f       = g * g
        denom_g = 1.0 - f
        denom_o = 1.0 - omega * f
        scaled_asymmetry[i] = (g - f) / denom_g if denom_g > 0.0 else 0.0
        scaled_ssa[i]       = (1.0 - f) * omega / denom_o if denom_o > 0.0 else 0.0
        scaled_tau[i]       = denom_o * tau

    # Cumulative and slant optical depths
    cumul_tau        = [0.0] * n_levels
    slant_tau        = [0.0] * n_levels
    effective_cosine = [1.0 / math.sqrt(LARGEST)] * n_levels

    if mu < 0.0:
        slant_tau[0] = _slant_optical_depth(0, nid[0], dsdh[0], scaled_tau)

    # Per-layer intermediate quantities
    lam             = [0.0] * n_layers
    gamma_ratio_arr = [0.0] * n_layers
    diffuse_cosine  = [0.0] * n_layers
    e1 = [0.0] * n_layers
    e2 = [0.0] * n_layers
    e3 = [0.0] * n_layers
    e4 = [0.0] * n_layers
    c_plus_top    = [0.0] * n_layers
    c_minus_top   = [0.0] * n_layers
    c_plus_bottom  = [0.0] * n_layers
    c_minus_bottom = [0.0] * n_layers

    for i in range(n_layers):
        g     = scaled_asymmetry[i]
        omega = scaled_ssa[i]
        # Clamp to avoid singularities (Toon et al. eq. 16)
        g     = math.copysign(min(abs(g), 1.0 - PRECISION), g)
        omega = min(omega, 1.0 - PRECISION)

        cumul_tau[i + 1] = cumul_tau[i] + scaled_tau[i]
        slant_tau[i + 1] = _slant_optical_depth(i + 1, nid[i + 1], dsdh[i + 1], scaled_tau)

        if nid[i + 1] >= 0:
            delta_slant = slant_tau[i + 1] - slant_tau[i]
            if delta_slant == 0.0:
                effective_cosine[i + 1] = math.sqrt(LARGEST)
            else:
                ec = (cumul_tau[i + 1] - cumul_tau[i]) / delta_slant
                effective_cosine[i + 1] = math.copysign(
                    max(abs(ec), 1.0 / math.sqrt(LARGEST)), ec)

        # Eddington gamma coefficients (Toon et al. 1989, Table 1, row 2)
        gamma1 =  (7.0 - omega * (4.0 + 3.0 * g)) / 4.0
        gamma2 = -(1.0 - omega * (4.0 - 3.0 * g)) / 4.0
        gamma3 =  (2.0 - 3.0 * g * mu) / 4.0
        gamma4 = 1.0 - gamma3
        diffuse_cosine[i] = 0.5

        lam[i] = math.sqrt(gamma1**2 - gamma2**2)
        gamma_ratio_arr[i] = (gamma1 - lam[i]) / gamma2 if gamma2 != 0.0 else 0.0

        exp_decay = math.exp(-lam[i] * scaled_tau[i])
        e1[i] = 1.0 + gamma_ratio_arr[i] * exp_decay
        e2[i] = 1.0 - gamma_ratio_arr[i] * exp_decay
        e3[i] = gamma_ratio_arr[i] + exp_decay
        e4[i] = gamma_ratio_arr[i] - exp_decay

        # Solar source functions (Toon et al. eqs. 23-24)
        beam_top    = math.exp(-slant_tau[i])
        beam_bottom = math.exp(-slant_tau[i + 1])
        divisor = lam[i]**2 - 1.0 / effective_cosine[i + 1]**2
        divisor = math.copysign(max(EPSILON, abs(divisor)), divisor)

        up_amp = omega * ((gamma1 - 1.0 / effective_cosine[i + 1]) * gamma3 + gamma4 * gamma2) / divisor
        dn_amp = omega * ((gamma1 + 1.0 / effective_cosine[i + 1]) * gamma4 + gamma2 * gamma3) / divisor

        c_plus_top[i]    = up_amp * beam_top
        c_minus_top[i]   = dn_amp * beam_top
        c_plus_bottom[i]  = up_amp * beam_bottom
        c_minus_bottom[i] = dn_amp * beam_bottom

    # Assemble tridiagonal system (Toon et al. eqs. 39-43)
    system_size = 2 * n_layers
    surface_solar_source = surface_albedo * mu * math.exp(-slant_tau[n_layers])

    lower_diag = [0.0] * (system_size - 1)
    main_diag  = [0.0] * system_size
    upper_diag = [0.0] * (system_size - 1)
    rhs        = [0.0] * system_size

    # Top boundary (Toon et al. eq. 39)
    main_diag[0]  = e1[0]
    upper_diag[0] = -e2[0]
    rhs[0]        = -c_minus_top[0]

    # Interior rows
    for k in range(n_layers - 1):
        even_row = 2 * k + 1
        lower_diag[even_row - 1] = e2[k + 1] * e1[k] - e3[k] * e4[k + 1]
        main_diag[even_row]      = e2[k] * e2[k + 1] - e4[k] * e4[k + 1]
        upper_diag[even_row]     = e1[k + 1] * e4[k + 1] - e2[k + 1] * e3[k + 1]
        rhs[even_row] = ((c_plus_top[k + 1] - c_plus_bottom[k]) * e2[k + 1]
                         - (c_minus_top[k + 1] - c_minus_bottom[k]) * e4[k + 1])

        odd_row = 2 * k + 2
        lower_diag[odd_row - 1] = e2[k] * e3[k] - e4[k] * e1[k]
        main_diag[odd_row]      = e1[k] * e1[k + 1] - e3[k] * e3[k + 1]
        upper_diag[odd_row]     = e3[k] * e4[k + 1] - e1[k] * e2[k + 1]
        rhs[odd_row] = (e3[k] * (c_plus_top[k + 1] - c_plus_bottom[k])
                        + e1[k] * (c_minus_bottom[k] - c_minus_top[k + 1]))

    # Bottom boundary (Toon et al. eq. 43)
    lower_diag[system_size - 2] = e1[n_layers - 1] - surface_albedo * e3[n_layers - 1]
    main_diag[system_size - 1]  = e2[n_layers - 1] - surface_albedo * e4[n_layers - 1]
    rhs[system_size - 1]        = (surface_solar_source
                                   - c_plus_bottom[n_layers - 1]
                                   + surface_albedo * c_minus_bottom[n_layers - 1])

    solution = _thomas_solve(lower_diag, main_diag, upper_diag, rhs)

    # Back-substitution (internal top-to-bottom order)
    direct_irradiance      = [0.0] * n_levels
    upwelling_irradiance   = [0.0] * n_levels
    downwelling_irradiance = [0.0] * n_levels
    direct_actinic_flux      = [0.0] * n_levels
    upwelling_actinic_flux   = [0.0] * n_levels
    downwelling_actinic_flux = [0.0] * n_levels

    direct_actinic_flux[0]      = math.exp(-slant_tau[0])
    direct_irradiance[0]        = mu * direct_actinic_flux[0]
    downwelling_irradiance[0]   = 0.0
    upwelling_irradiance[0]     = solution[0] * e3[0] - solution[1] * e4[0] + c_plus_top[0]
    downwelling_actinic_flux[0] = downwelling_irradiance[0] / diffuse_cosine[0]
    upwelling_actinic_flux[0]   = upwelling_irradiance[0] / diffuse_cosine[0]

    for lev in range(1, n_levels):
        layer   = lev - 1
        row_idx = 2 * layer
        direct_actinic_flux[lev]      = math.exp(-slant_tau[lev])
        direct_irradiance[lev]        = mu * direct_actinic_flux[lev]
        downwelling_irradiance[lev]   = (solution[row_idx] * e3[layer]
                                         + solution[row_idx + 1] * e4[layer]
                                         + c_minus_bottom[layer])
        upwelling_irradiance[lev]     = (solution[row_idx] * e1[layer]
                                         + solution[row_idx + 1] * e2[layer]
                                         + c_plus_bottom[layer])
        downwelling_actinic_flux[lev] = downwelling_irradiance[lev] / diffuse_cosine[layer]
        upwelling_actinic_flux[lev]   = upwelling_irradiance[lev] / diffuse_cosine[layer]

    # Reverse internal top-to-bottom → output bottom-to-top
    return {
        "direct_irradiance":        list(reversed(direct_irradiance)),
        "upwelling_irradiance":     list(reversed(upwelling_irradiance)),
        "downwelling_irradiance":   list(reversed(downwelling_irradiance)),
        "direct_actinic_flux":      list(reversed(direct_actinic_flux)),
        "upwelling_actinic_flux":   list(reversed(upwelling_actinic_flux)),
        "downwelling_actinic_flux": list(reversed(downwelling_actinic_flux)),
    }


# ── CSV output ─────────────────────────────────────────────────────────────────

FIELDS = [
    "level",
    "direct_irradiance",
    "upwelling_irradiance",
    "downwelling_irradiance",
    "direct_actinic_flux",
    "upwelling_actinic_flux",
    "downwelling_actinic_flux",
]

REFERENCE_DIR = os.path.join(os.path.dirname(__file__), "reference")


def write_csv(sza_rad, result):
    filename = os.path.join(REFERENCE_DIR, f"delta_eddington_sza_{sza_rad:.4f}.csv")
    n_levels = len(result["direct_irradiance"])
    with open(filename, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(
            ["# atmosphere: 3-layer, tau=0.5/layer, ssa=0.9/layer, g=0.85/layer, albedo=0.1"]
        )
        writer.writerow([f"# sza_rad={sza_rad}"])
        writer.writerow(["# level 0 = ground, level n_layers = TOA"])
        writer.writerow(FIELDS)
        for lev in range(n_levels):
            writer.writerow([
                lev,
                f"{result['direct_irradiance'][lev]:.15e}",
                f"{result['upwelling_irradiance'][lev]:.15e}",
                f"{result['downwelling_irradiance'][lev]:.15e}",
                f"{result['direct_actinic_flux'][lev]:.15e}",
                f"{result['upwelling_actinic_flux'][lev]:.15e}",
                f"{result['downwelling_actinic_flux'][lev]:.15e}",
            ])
    print(f"Wrote {filename}")


if __name__ == "__main__":
    os.makedirs(REFERENCE_DIR, exist_ok=True)
    for sza in SZAS:
        result = solve(
            sza_rad=sza,
            optical_depth=OPTICAL_DEPTH,
            single_scattering_albedo=SINGLE_SCATTERING_ALBEDO,
            asymmetry_parameter=ASYMMETRY_PARAMETER,
            surface_albedo=SURFACE_ALBEDO,
            alt_edges_m=ALT_EDGES_M,
        )
        write_csv(sza, result)
    print("Done.")
