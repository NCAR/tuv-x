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
#include <ranges>
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
    constexpr double largest   = 1.0e36;
    constexpr double precision = 1.0e-7;
    constexpr double epsilon   = 1.0e-3;

    const std::size_t n_columns = solar_zenith_angles.size();

    const auto& vertical_grid   = grids.at("altitude [m]");
    const auto& wavelength_grid = grids.at("wavelength [m]");
    const auto& surface_albedo  = profiles.at("surface albedo [1]");

    const std::size_t n_layers      = vertical_grid.NumberOfSections();
    const std::size_t n_levels      = n_layers + 1;
    const std::size_t n_wavelengths = wavelength_grid.NumberOfSections();
    const std::size_t system_size   = 2 * n_layers;

    for (std::size_t col = 0; col < n_columns; ++col)
    {
      const double cosine_solar_zenith = std::cos(static_cast<double>(solar_zenith_angles[col]));

      const std::size_t grid_col = (vertical_grid.NumberOfColumns() == 1) ? 0 : col;
      std::vector<double> altitude_edges(n_levels);
      std::ranges::transform(
          std::views::iota(std::size_t{ 0 }, n_levels),
          altitude_edges.begin(),
          [&](std::size_t lev) { return vertical_grid.edges_(lev, grid_col); });

      SphericalGeometry geometry;
      geometry.SetParameters(solar_zenith_angles[col], altitude_edges);

      const std::size_t albedo_column = (surface_albedo.mid_point_values_.Size2() == 1) ? 0 : col;

      for (std::size_t wl = 0; wl < n_wavelengths; ++wl)
      {
        const double surface_reflectivity = surface_albedo.mid_point_values_(wl, albedo_column);

        // ── Delta-scaling (Toon et al. 1989, Section 3) ───────────────────
        std::vector<double> scaled_asymmetry(n_layers);
        std::vector<double> scaled_single_scattering_albedo(n_layers);
        std::vector<double> scaled_optical_depth(n_layers);
        for (std::size_t i = 0; i < n_layers; ++i)
        {
          const double g     = constituent_states.asymmetry_parameter_(wl, i, col);
          const double omega = constituent_states.single_scattering_albedo_(wl, i, col);
          const double tau   = constituent_states.optical_depth_(wl, i, col);
          const double f            = g * g;
          const double denom_g      = 1.0 - f;
          const double denom_omega  = 1.0 - omega * f;
          scaled_asymmetry[i]                = (denom_g > 0.0) ? (g - f) / denom_g : 0.0;
          scaled_single_scattering_albedo[i] = (denom_omega > 0.0) ? (1.0 - f) * omega / denom_omega : 0.0;
          scaled_optical_depth[i]            = denom_omega * tau;
        }

        // ── Cumulative and slant optical depths ───────────────────────────
        std::vector<double> cumulative_optical_depth(n_levels, 0.0);
        std::vector<double> slant_optical_depths(n_levels, 0.0);
        // effective_cosine[level]: effective cosine of the direct beam in the
        // layer below that level; initialised to near-zero (large slant path).
        std::vector<double> effective_cosine(n_levels, 1.0 / std::sqrt(largest));

        // For a below-horizon sun, the slant path at TOA is non-zero because the
        // beam arrives from below; compute the initial slant optical depth there.
        if (cosine_solar_zenith < 0.0)
        {
          slant_optical_depths[0] = SlantOpticalDepth(
              0, geometry.nid_[0], geometry.dsdh_[0], scaled_optical_depth);
        }

        // Per-layer intermediate arrays
        std::vector<double> lambda(n_layers), gamma_ratio(n_layers), diffuse_cosine(n_layers);
        std::vector<double> e1(n_layers), e2(n_layers), e3(n_layers), e4(n_layers);
        std::vector<double> c_plus_top(n_layers), c_minus_top(n_layers);
        std::vector<double> c_plus_bottom(n_layers), c_minus_bottom(n_layers);

        for (std::size_t i = 0; i < n_layers; ++i)
        {
          // Clamp to avoid singularities (Toon et al. eq. 16)
          double g     = scaled_asymmetry[i];
          double omega = scaled_single_scattering_albedo[i];
          g     = std::copysign(std::min(std::abs(g), 1.0 - precision), g);
          omega = std::min(omega, 1.0 - precision);

          cumulative_optical_depth[i + 1] = cumulative_optical_depth[i] + scaled_optical_depth[i];

          slant_optical_depths[i + 1] = SlantOpticalDepth(
              i + 1, geometry.nid_[i + 1], geometry.dsdh_[i + 1], scaled_optical_depth);

          if (geometry.nid_[i + 1].has_value())
          {
            const double delta_slant = slant_optical_depths[i + 1] - slant_optical_depths[i];
            if (delta_slant == 0.0)
            {
              effective_cosine[i + 1] = std::sqrt(largest);
            }
            else
            {
              effective_cosine[i + 1] =
                  (cumulative_optical_depth[i + 1] - cumulative_optical_depth[i]) / delta_slant;
              effective_cosine[i + 1] = std::copysign(
                  std::max(std::abs(effective_cosine[i + 1]), 1.0 / std::sqrt(largest)),
                  effective_cosine[i + 1]);
            }
          }

          // Eddington gamma coefficients (Toon et al. 1989, Table 1, row 2)
          const double gamma1 =  (7.0 - omega * (4.0 + 3.0 * g)) / 4.0;
          const double gamma2 = -(1.0 - omega * (4.0 - 3.0 * g)) / 4.0;
          const double gamma3 =  (2.0 - 3.0 * g * cosine_solar_zenith) / 4.0;
          const double gamma4 = 1.0 - gamma3;
          diffuse_cosine[i] = 0.5;  // Eddington approximation: <cos theta> = 1/2

          lambda[i]     = std::sqrt(gamma1 * gamma1 - gamma2 * gamma2);
          gamma_ratio[i] = (gamma2 != 0.0) ? (gamma1 - lambda[i]) / gamma2 : 0.0;

          const double exp_decay = std::exp(-lambda[i] * scaled_optical_depth[i]);
          e1[i] = 1.0 + gamma_ratio[i] * exp_decay;
          e2[i] = 1.0 - gamma_ratio[i] * exp_decay;
          e3[i] = gamma_ratio[i] + exp_decay;
          e4[i] = gamma_ratio[i] - exp_decay;

          // Solar source functions C+ and C- (Toon et al. eqs. 23–24)
          const double beam_at_top    = std::exp(-slant_optical_depths[i]);
          const double beam_at_bottom = std::exp(-slant_optical_depths[i + 1]);

          double divisor = lambda[i] * lambda[i] - 1.0 / (effective_cosine[i + 1] * effective_cosine[i + 1]);
          divisor = std::copysign(std::max(epsilon, std::abs(divisor)), divisor);

          const double upward_amplitude   = omega * ((gamma1 - 1.0 / effective_cosine[i + 1]) * gamma3 + gamma4 * gamma2) / divisor;
          const double downward_amplitude = omega * ((gamma1 + 1.0 / effective_cosine[i + 1]) * gamma4 + gamma2 * gamma3) / divisor;

          c_plus_top[i]    = upward_amplitude * beam_at_top;
          c_minus_top[i]   = downward_amplitude * beam_at_top;
          c_plus_bottom[i] = upward_amplitude * beam_at_bottom;
          c_minus_bottom[i] = downward_amplitude * beam_at_bottom;
        }

        // ── Assemble tridiagonal system (Toon et al. eqs. 39–43) ──────────
        const double surface_solar_source =
            surface_reflectivity * cosine_solar_zenith * std::exp(-slant_optical_depths[n_layers]);

        TridiagonalMatrix<double> tridiagonal(system_size);
        Array1D<double>           rhs(system_size);

        // Top boundary: zero diffuse incidence (Toon et al. eq. 39)
        tridiagonal.main_diagonal_[0]  = e1[0];
        tridiagonal.upper_diagonal_[0] = -e2[0];
        rhs[0]                         = -c_minus_top[0];

        // Interior rows: flux continuity at layer interfaces
        for (std::size_t k = 0; k < n_layers - 1; ++k)
        {
          const std::size_t even_row = 2 * k + 1;
          tridiagonal.lower_diagonal_[even_row - 1] = e2[k + 1] * e1[k] - e3[k] * e4[k + 1];
          tridiagonal.main_diagonal_[even_row]      = e2[k] * e2[k + 1] - e4[k] * e4[k + 1];
          tridiagonal.upper_diagonal_[even_row]     = e1[k + 1] * e4[k + 1] - e2[k + 1] * e3[k + 1];
          rhs[even_row] = (c_plus_top[k + 1] - c_plus_bottom[k]) * e2[k + 1]
                        - (c_minus_top[k + 1] - c_minus_bottom[k]) * e4[k + 1];

          const std::size_t odd_row = 2 * k + 2;
          tridiagonal.lower_diagonal_[odd_row - 1] = e2[k] * e3[k] - e4[k] * e1[k];
          tridiagonal.main_diagonal_[odd_row]      = e1[k] * e1[k + 1] - e3[k] * e3[k + 1];
          tridiagonal.upper_diagonal_[odd_row]     = e3[k] * e4[k + 1] - e1[k] * e2[k + 1];
          rhs[odd_row] = e3[k] * (c_plus_top[k + 1] - c_plus_bottom[k])
                       + e1[k] * (c_minus_bottom[k] - c_minus_top[k + 1]);
        }

        // Bottom boundary: surface albedo condition (Toon et al. eq. 43)
        tridiagonal.lower_diagonal_[system_size - 2] =
            e1[n_layers - 1] - surface_reflectivity * e3[n_layers - 1];
        tridiagonal.main_diagonal_[system_size - 1] =
            e2[n_layers - 1] - surface_reflectivity * e4[n_layers - 1];
        rhs[system_size - 1] =
            surface_solar_source - c_plus_bottom[n_layers - 1] + surface_reflectivity * c_minus_bottom[n_layers - 1];

        tuvx::Solve(tridiagonal, rhs);
        const Array1D<double>& solution = rhs;

        // ── Back-substitute: compute fluxes (internal top-to-bottom order) ─
        std::vector<double> direct_irradiance(n_levels);
        std::vector<double> upwelling_irradiance(n_levels);
        std::vector<double> downwelling_irradiance(n_levels);
        std::vector<double> direct_actinic_flux(n_levels);
        std::vector<double> upwelling_actinic_flux(n_levels);
        std::vector<double> downwelling_actinic_flux(n_levels);

        // Level 0 = TOA
        direct_actinic_flux[0]      = std::exp(-slant_optical_depths[0]);
        direct_irradiance[0]        = cosine_solar_zenith * direct_actinic_flux[0];
        downwelling_irradiance[0]   = 0.0;
        upwelling_irradiance[0]     = solution[0] * e3[0] - solution[1] * e4[0] + c_plus_top[0];
        downwelling_actinic_flux[0] = downwelling_irradiance[0] / diffuse_cosine[0];
        upwelling_actinic_flux[0]   = upwelling_irradiance[0] / diffuse_cosine[0];

        for (std::size_t lev = 1; lev <= n_layers; ++lev)
        {
          const std::size_t layer   = lev - 1;
          const std::size_t row_idx = 2 * layer;
          direct_actinic_flux[lev]      = std::exp(-slant_optical_depths[lev]);
          direct_irradiance[lev]        = cosine_solar_zenith * direct_actinic_flux[lev];
          downwelling_irradiance[lev]   = solution[row_idx] * e3[layer] + solution[row_idx + 1] * e4[layer] + c_minus_bottom[layer];
          upwelling_irradiance[lev]     = solution[row_idx] * e1[layer] + solution[row_idx + 1] * e2[layer] + c_plus_bottom[layer];
          downwelling_actinic_flux[lev] = downwelling_irradiance[lev] / diffuse_cosine[layer];
          upwelling_actinic_flux[lev]   = upwelling_irradiance[lev] / diffuse_cosine[layer];
        }

        // ── Write results reversed to bottom-to-top convention ─────────────
        for (std::size_t lev = 0; lev <= n_layers; ++lev)
        {
          const std::size_t out_level = n_layers - lev;
          radiation_field.spectral_irradiance_.direct_(wl, out_level, col)      = direct_irradiance[lev];
          radiation_field.spectral_irradiance_.upwelling_(wl, out_level, col)   = upwelling_irradiance[lev];
          radiation_field.spectral_irradiance_.downwelling_(wl, out_level, col) = downwelling_irradiance[lev];
          radiation_field.actinic_flux_.direct_(wl, out_level, col)             = direct_actinic_flux[lev];
          radiation_field.actinic_flux_.upwelling_(wl, out_level, col)          = upwelling_actinic_flux[lev];
          radiation_field.actinic_flux_.downwelling_(wl, out_level, col)        = downwelling_actinic_flux[lev];
        }
      }  // wavelength loop
    }  // column loop
  }

}  // namespace tuvx
