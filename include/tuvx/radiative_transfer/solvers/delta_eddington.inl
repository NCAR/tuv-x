// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
//
// Ported from src/radiative_transfer/solvers/delta_eddington.F90
// Reference: Toon et al., J. Geophys. Res. 94(D13), 16287–16301, 1989.
//            doi:10.1029/JD094iD13p16287
#include <tuvx/linear_algebra/linear_algebra.hpp>
#include <tuvx/radiative_transfer/spherical_geometry.hpp>

#include <algorithm>
#include <cmath>
#include <vector>

namespace tuvx
{

  // ──────────────────────────────────────────────────────────────────────────
  // Index conventions used throughout this file:
  //   Layer index 0 = topmost layer (just below TOA); n_layers-1 = bottommost.
  //   Level index 0 = TOA; n_layers = ground.
  //   Altitude grid edges are stored bottom-to-top (index 0 = ground).
  //   Radiation field arrays are stored bottom-to-top (level 0 = ground).
  // ──────────────────────────────────────────────────────────────────────────

  template<
      typename T,
      typename GridPolicy,
      typename ProfilePolicy,
      typename ConstituentStatePolicy,
      typename RadiationFieldPolicy>
  inline void DeltaEddington::Solve(
      const std::vector<T>&                       solar_zenith_angles,
      const std::map<std::string, GridPolicy>&    grids,
      const std::map<std::string, ProfilePolicy>& profiles,
      const ConstituentStatePolicy&               constituent_states,
      RadiationFieldPolicy&                       radiation_field) const
  {
    constexpr double kLargest = 1.0e36;
    constexpr double kPrecis  = 1.0e-7;
    constexpr double kEps     = 1.0e-3;

    const std::size_t n_columns = solar_zenith_angles.size();

    const auto& vertical_grid   = grids.at("altitude [m]");
    const auto& wavelength_grid = grids.at("wavelength [m]");
    const auto& surface_albedo  = profiles.at("surface albedo [1]");

    const std::size_t n_layers      = vertical_grid.NumberOfSections();
    const std::size_t n_levels      = n_layers + 1;
    const std::size_t n_wavelengths = wavelength_grid.NumberOfSections();
    const std::size_t mrows         = 2 * n_layers;

    for (std::size_t col = 0; col < n_columns; ++col)
    {
      const double mu = std::cos(static_cast<double>(solar_zenith_angles[col]));

      // Collect altitude edges for this column (bottom-to-top)
      const std::size_t grid_col = (vertical_grid.NumberOfColumns() == 1) ? 0 : col;
      std::vector<double> altitude_edges(n_levels);
      for (std::size_t lev = 0; lev < n_levels; ++lev)
      {
        altitude_edges[lev] = vertical_grid.edges_(lev, grid_col);
      }

      SphericalGeometry sph_geom;
      sph_geom.SetParameters(static_cast<double>(solar_zenith_angles[col]), altitude_edges);

      const std::size_t alb_col = (surface_albedo.mid_point_values_.Size2() == 1) ? 0 : col;

      for (std::size_t wl = 0; wl < n_wavelengths; ++wl)
      {
        const double rsfc = surface_albedo.mid_point_values_(wl, alb_col);

        // Unpack optical properties for this column+wavelength (top-to-bottom)
        std::vector<double> tauu(n_layers), omu(n_layers), gu(n_layers);
        for (std::size_t layer = 0; layer < n_layers; ++layer)
        {
          tauu[layer] = constituent_states.optical_depth_(wl, layer, col);
          omu[layer]  = constituent_states.single_scattering_albedo_(wl, layer, col);
          gu[layer]   = constituent_states.asymmetry_parameter_(wl, layer, col);
        }

        // ── Delta-scaling ─────────────────────────────────────────────────
        std::vector<double> gi(n_layers), omi(n_layers), taun(n_layers);
        for (std::size_t i = 0; i < n_layers; ++i)
        {
          const double f = gu[i] * gu[i];
          const double denom_g = 1.0 - f;
          gi[i]   = (denom_g > 0.0) ? (gu[i] - f) / denom_g : 0.0;
          const double denom_o = 1.0 - omu[i] * f;
          omi[i]  = (denom_o > 0.0) ? (1.0 - f) * omu[i] / denom_o : 0.0;
          taun[i] = (1.0 - omu[i] * f) * tauu[i];
        }

        // ── Cumulative and slant optical depths ───────────────────────────
        std::vector<double> tauc(n_levels, 0.0);
        std::vector<double> tausla(n_levels, 0.0);
        // mu2[level]: effective cosine of the direct beam in the layer below
        // that level; initialised to near-zero (large slant path).
        std::vector<double> mu2(n_levels, 1.0 / std::sqrt(kLargest));

        if (mu < 0.0)
        {
          tausla[0] = SlantOpticalDepth(0, sph_geom.nid_[0], sph_geom.dsdh_[0], taun);
        }

        // Per-layer intermediate arrays
        std::vector<double> lam(n_layers), bgam(n_layers), mu1(n_layers);
        std::vector<double> e1(n_layers), e2(n_layers), e3(n_layers), e4(n_layers);
        std::vector<double> cup(n_layers), cdn(n_layers);
        std::vector<double> cuptn(n_layers), cdntn(n_layers);

        for (std::size_t i = 0; i < n_layers; ++i)
        {
          // Clamp to avoid singularities (Toon et al. eq. 16)
          double g  = gi[i];
          double om = omi[i];
          const double tempg = std::min(std::abs(g), 1.0 - kPrecis);
          g  = std::copysign(tempg, g);
          om = std::min(om, 1.0 - kPrecis);

          tauc[i + 1] = tauc[i] + taun[i];

          tausla[i + 1] = SlantOpticalDepth(
              i + 1, sph_geom.nid_[i + 1], sph_geom.dsdh_[i + 1], taun);

          if (sph_geom.nid_[i + 1] >= 0)
          {
            const double delta_tausla = tausla[i + 1] - tausla[i];
            if (delta_tausla == 0.0)
            {
              mu2[i + 1] = std::sqrt(kLargest);
            }
            else
            {
              mu2[i + 1] = (tauc[i + 1] - tauc[i]) / delta_tausla;
              mu2[i + 1] = std::copysign(
                  std::max(std::abs(mu2[i + 1]), 1.0 / std::sqrt(kLargest)), mu2[i + 1]);
            }
          }

          // Eddington gamma coefficients (Toon et al. 1989, Table 1, row 2)
          const double gam1 =  (7.0 - om * (4.0 + 3.0 * g)) / 4.0;
          const double gam2 = -(1.0 - om * (4.0 - 3.0 * g)) / 4.0;
          const double gam3 =  (2.0 - 3.0 * g * mu) / 4.0;
          const double gam4 = 1.0 - gam3;
          mu1[i] = 0.5;

          lam[i]  = std::sqrt(gam1 * gam1 - gam2 * gam2);
          bgam[i] = (gam2 != 0.0) ? (gam1 - lam[i]) / gam2 : 0.0;

          const double expon = std::exp(-lam[i] * taun[i]);
          e1[i] = 1.0 + bgam[i] * expon;
          e2[i] = 1.0 - bgam[i] * expon;
          e3[i] = bgam[i] + expon;
          e4[i] = bgam[i] - expon;

          // Solar source functions C+ and C- (Toon et al. eqs. 23–24)
          const double expon0 = std::exp(-tausla[i]);
          const double expon1 = std::exp(-tausla[i + 1]);

          double divisr = lam[i] * lam[i] - 1.0 / (mu2[i + 1] * mu2[i + 1]);
          divisr = std::copysign(std::max(kEps, std::abs(divisr)), divisr);

          const double up = om * ((gam1 - 1.0 / mu2[i + 1]) * gam3 + gam4 * gam2) / divisr;
          const double dn = om * ((gam1 + 1.0 / mu2[i + 1]) * gam4 + gam2 * gam3) / divisr;

          cup[i]   = up * expon0;
          cdn[i]   = dn * expon0;
          cuptn[i] = up * expon1;
          cdntn[i] = dn * expon1;
        }

        // ── Assemble tridiagonal system (Toon et al. eqs. 39–43) ──────────
        const double ssfc = rsfc * mu * std::exp(-tausla[n_layers]);

        TridiagonalMatrix<double> tri(mrows);
        Array1D<double>           rhs(mrows);

        // Top boundary: zero diffuse incidence (Toon et al. eq. 39)
        tri.main_diagonal_[0]  = e1[0];
        tri.upper_diagonal_[0] = -e2[0];
        rhs[0]                 = -cdn[0];  // fdn0 = 0

        // Interior rows: flux continuity at layer interfaces
        for (std::size_t k = 0; k < n_layers - 1; ++k)
        {
          const std::size_t row_e = 2 * k + 1;
          tri.lower_diagonal_[row_e - 1] = e2[k + 1] * e1[k] - e3[k] * e4[k + 1];
          tri.main_diagonal_[row_e]      = e2[k] * e2[k + 1] - e4[k] * e4[k + 1];
          tri.upper_diagonal_[row_e]     = e1[k + 1] * e4[k + 1] - e2[k + 1] * e3[k + 1];
          rhs[row_e] = (cup[k + 1] - cuptn[k]) * e2[k + 1]
                     - (cdn[k + 1] - cdntn[k]) * e4[k + 1];

          const std::size_t row_o = 2 * k + 2;
          tri.lower_diagonal_[row_o - 1] = e2[k] * e3[k] - e4[k] * e1[k];
          tri.main_diagonal_[row_o]      = e1[k] * e1[k + 1] - e3[k] * e3[k + 1];
          tri.upper_diagonal_[row_o]     = e3[k] * e4[k + 1] - e1[k] * e2[k + 1];
          rhs[row_o] = e3[k] * (cup[k + 1] - cuptn[k])
                     + e1[k] * (cdntn[k] - cdn[k + 1]);
        }

        // Bottom boundary: surface albedo condition (Toon et al. eq. 43)
        tri.lower_diagonal_[mrows - 2] = e1[n_layers - 1] - rsfc * e3[n_layers - 1];
        tri.main_diagonal_[mrows - 1]  = e2[n_layers - 1] - rsfc * e4[n_layers - 1];
        rhs[mrows - 1] = ssfc - cuptn[n_layers - 1] + rsfc * cdntn[n_layers - 1];

        tuvx::Solve(tri, rhs);  // solution overwrites rhs (renamed y below for clarity)
        const Array1D<double>& y = rhs;

        // ── Back-substitute: compute fluxes (internal top-to-bottom order) ─
        std::vector<double> edr(n_levels), eup_v(n_levels), edn_v(n_levels);
        std::vector<double> fdr(n_levels), fup_v(n_levels), fdn_v(n_levels);

        // Level 0 = TOA
        fdr[0]   = std::exp(-tausla[0]);
        edr[0]   = mu * fdr[0];
        edn_v[0] = 0.0;
        eup_v[0] = y[0] * e3[0] - y[1] * e4[0] + cup[0];
        fdn_v[0] = edn_v[0] / mu1[0];
        fup_v[0] = eup_v[0] / mu1[0];

        for (std::size_t lev = 1; lev <= n_layers; ++lev)
        {
          const std::size_t j       = lev - 1;  // 0-based layer
          const std::size_t row_idx = 2 * j;    // index into y
          fdr[lev]   = std::exp(-tausla[lev]);
          edr[lev]   = mu * fdr[lev];
          edn_v[lev] = y[row_idx] * e3[j] + y[row_idx + 1] * e4[j] + cdntn[j];
          eup_v[lev] = y[row_idx] * e1[j] + y[row_idx + 1] * e2[j] + cuptn[j];
          fdn_v[lev] = edn_v[lev] / mu1[j];
          fup_v[lev] = eup_v[lev] / mu1[j];
        }

        // ── Write results reversed to bottom-to-top convention ─────────────
        for (std::size_t lev = 0; lev <= n_layers; ++lev)
        {
          const std::size_t out = n_layers - lev;  // flip: internal 0→out n_layers
          radiation_field.spectral_irradiance_.direct_(wl, out, col)      = edr[lev];
          radiation_field.spectral_irradiance_.upwelling_(wl, out, col)   = eup_v[lev];
          radiation_field.spectral_irradiance_.downwelling_(wl, out, col) = edn_v[lev];
          radiation_field.actinic_flux_.direct_(wl, out, col)             = fdr[lev];
          radiation_field.actinic_flux_.upwelling_(wl, out, col)          = fup_v[lev];
          radiation_field.actinic_flux_.downwelling_(wl, out, col)        = fdn_v[lev];
        }
      }  // wavelength loop
    }  // column loop
  }

}  // namespace tuvx
