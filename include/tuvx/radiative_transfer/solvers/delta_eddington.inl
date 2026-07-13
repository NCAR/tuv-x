// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
//
// Delta-Eddington two-stream radiative transfer solver.
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

  // ── Helper: per-layer result bundle ──────────────────────────────────────

  struct LayerCoefficients
  {
    double lambda_;
    double gamma_ratio_;
    double diffuse_cosine_;
    double e1_;
    double e2_;
    double e3_;
    double e4_;
    double c_plus_top_;
    double c_minus_top_;
    double c_plus_bottom_;
    double c_minus_bottom_;
  };

  // ── Helper: effective cosine of the direct beam in a layer ───────────────
  // Computes the effective cosine without deeply nested conditionals inside Solve.
  inline double ComputeEffectiveCosine(bool nid_has_value, double delta_od, double delta_slant, double large_sentinel)
  {
    const double INV_SQRT = 1.0 / std::sqrt(large_sentinel);
    if (!nid_has_value)
    {
      return INV_SQRT;
    }
    if (delta_slant == 0.0)
    {
      return std::sqrt(large_sentinel);
    }
    const double EFF = delta_od / delta_slant;
    return std::copysign(std::max(std::abs(EFF), INV_SQRT), EFF);
  }

  // ── Helper: per-layer Eddington coefficients and source terms ─────────────
  inline LayerCoefficients ComputeLayerCoefficients(
      double g_scaled,
      double omega_scaled,
      double tau_scaled,
      double cosine_solar_zenith,
      double effective_cosine_layer,
      double beam_at_top,
      double beam_at_bottom,
      double singularity_guard,
      double source_divisor_floor)
  {
    // Clamp to avoid singularities (Toon et al. eq. 16)
    double g = std::copysign(std::min(std::abs(g_scaled), 1.0 - singularity_guard), g_scaled);
    double omega = std::min(omega_scaled, 1.0 - singularity_guard);

    // Eddington gamma coefficients (Toon et al. 1989, Table 1, row 2)
    const double GAMMA1 = (7.0 - (omega * (4.0 + (3.0 * g)))) / 4.0;
    const double GAMMA2 = -(1.0 - (omega * (4.0 - (3.0 * g)))) / 4.0;
    const double GAMMA3 = (2.0 - (3.0 * g * cosine_solar_zenith)) / 4.0;
    const double GAMMA4 = 1.0 - GAMMA3;

    const double LAM = std::sqrt((GAMMA1 * GAMMA1) - (GAMMA2 * GAMMA2));
    const double GAMMA_RATIO = (GAMMA2 != 0.0) ? (GAMMA1 - LAM) / GAMMA2 : 0.0;
    const double EXP_DECAY = std::exp(-(LAM * tau_scaled));

    // Solar source functions C+ and C- (Toon et al. eqs. 23–24)
    double divisor = (LAM * LAM) - (1.0 / (effective_cosine_layer * effective_cosine_layer));
    divisor = std::copysign(std::max(source_divisor_floor, std::abs(divisor)), divisor);

    const double UPWARD_AMP = omega * (((GAMMA1 - (1.0 / effective_cosine_layer)) * GAMMA3) + (GAMMA4 * GAMMA2)) / divisor;
    const double DOWNWARD_AMP = omega * (((GAMMA1 + (1.0 / effective_cosine_layer)) * GAMMA4) + (GAMMA2 * GAMMA3)) / divisor;

    return { .lambda_ = LAM,
             .gamma_ratio_ = GAMMA_RATIO,
             .diffuse_cosine_ = 0.5,  // Eddington approximation: <cos theta> = 1/2
             .e1_ = 1.0 + (GAMMA_RATIO * EXP_DECAY),
             .e2_ = 1.0 - (GAMMA_RATIO * EXP_DECAY),
             .e3_ = GAMMA_RATIO + EXP_DECAY,
             .e4_ = GAMMA_RATIO - EXP_DECAY,
             .c_plus_top_ = UPWARD_AMP * beam_at_top,
             .c_minus_top_ = DOWNWARD_AMP * beam_at_top,
             .c_plus_bottom_ = UPWARD_AMP * beam_at_bottom,
             .c_minus_bottom_ = DOWNWARD_AMP * beam_at_bottom };
  }

  template<
      typename T,
      typename GridPolicy,
      typename ProfilePolicy,
      typename ConstituentStatePolicy,
      typename RadiationFieldPolicy>
  inline void DeltaEddington::Solve(
      const std::vector<T>& solar_zenith_angles,
      const std::map<std::string, GridPolicy>& grids,
      const std::map<std::string, ProfilePolicy>& profiles,
      const ConstituentStatePolicy& constituent_states,
      RadiationFieldPolicy& radiation_field) const
  {
    constexpr double LARGE_COSINE_SENTINEL = 1.0e36;
    constexpr double SINGULARITY_GUARD = 1.0e-7;
    constexpr double SOURCE_DIVISOR_FLOOR = 1.0e-3;

    const std::size_t N_COLUMNS = solar_zenith_angles.size();

    const auto& vertical_grid = grids.at("altitude [m]");
    const auto& wavelength_grid = grids.at("wavelength [m]");
    const auto& surface_albedo = profiles.at("surface albedo [1]");

    const std::size_t N_LAYERS = vertical_grid.NumberOfSections();
    const std::size_t N_LEVELS = N_LAYERS + 1;
    const std::size_t N_WAVELENGTHS = wavelength_grid.NumberOfSections();
    const std::size_t SYSTEM_SIZE = 2 * N_LAYERS;

    for (std::size_t col = 0; col < N_COLUMNS; ++col)
    {
      const double COSINE_SOLAR_ZENITH = std::cos(static_cast<double>(solar_zenith_angles[col]));

      const std::size_t GRID_COL = (vertical_grid.NumberOfColumns() == 1) ? 0 : col;
      std::vector<double> altitude_edges(N_LEVELS);
      std::ranges::transform(
          std::views::iota(std::size_t{ 0 }, N_LEVELS),
          altitude_edges.begin(),
          [&](std::size_t lev) { return vertical_grid.edges_(lev, GRID_COL); });

      SphericalGeometry geom;
      geom.SetParameters(static_cast<double>(solar_zenith_angles[col]), altitude_edges);

      const std::size_t ALBEDO_COLUMN = (surface_albedo.mid_point_values_.Size2() == 1) ? 0 : col;

      for (std::size_t wl = 0; wl < N_WAVELENGTHS; ++wl)
      {
        const double SURFACE_REFLECTIVITY = surface_albedo.mid_point_values_(wl, ALBEDO_COLUMN);

        // ── Delta-scaling (Toon et al. 1989, Section 3) ───────────────────
        std::vector<double> scaled_asymmetry(N_LAYERS);
        std::vector<double> scaled_single_scattering_albedo(N_LAYERS);
        std::vector<double> scaled_optical_depth(N_LAYERS);
        for (std::size_t i = 0; i < N_LAYERS; ++i)
        {
          const double G = constituent_states.asymmetry_parameter_(wl, i, col);
          const double OMEGA = constituent_states.single_scattering_albedo_(wl, i, col);
          const double TAU = constituent_states.optical_depth_(wl, i, col);
          const double F = G * G;
          const double DENOM_G = 1.0 - F;
          const double DENOM_OMEGA = 1.0 - (OMEGA * F);
          scaled_asymmetry[i] = (DENOM_G > 0.0) ? (G - F) / DENOM_G : 0.0;
          scaled_single_scattering_albedo[i] = (DENOM_OMEGA > 0.0) ? ((1.0 - F) * OMEGA) / DENOM_OMEGA : 0.0;
          scaled_optical_depth[i] = DENOM_OMEGA * TAU;
        }

        // ── Cumulative and slant optical depths ───────────────────────────
        std::vector<double> cumulative_optical_depth(N_LEVELS, 0.0);
        std::vector<double> slant_optical_depths(N_LEVELS, 0.0);
        // effective_cosine[level]: effective cosine of the direct beam in the
        // layer below that level; initialised to near-zero (large slant path).
        std::vector<double> effective_cosine(N_LEVELS, 1.0 / std::sqrt(LARGE_COSINE_SENTINEL));

        // For a below-horizon sun, the slant path at TOA is non-zero because the
        // beam arrives from below; compute the initial slant optical depth there.
        if (COSINE_SOLAR_ZENITH < 0.0)
        {
          slant_optical_depths[0] = SlantOpticalDepth(0, geom.nid_[0], geom.dsdh_[0], scaled_optical_depth);
        }

        // Per-layer intermediate arrays
        std::vector<double> lambda(N_LAYERS), gamma_ratio(N_LAYERS), diffuse_cosine(N_LAYERS);
        std::vector<double> e1(N_LAYERS), e2(N_LAYERS), e3(N_LAYERS), e4(N_LAYERS);
        std::vector<double> c_plus_top(N_LAYERS), c_minus_top(N_LAYERS);
        std::vector<double> c_plus_bottom(N_LAYERS), c_minus_bottom(N_LAYERS);

        for (std::size_t i = 0; i < N_LAYERS; ++i)
        {
          cumulative_optical_depth[i + 1] = cumulative_optical_depth[i] + scaled_optical_depth[i];

          slant_optical_depths[i + 1] = SlantOpticalDepth(i + 1, geom.nid_[i + 1], geom.dsdh_[i + 1], scaled_optical_depth);

          effective_cosine[i + 1] = ComputeEffectiveCosine(
              geom.nid_[i + 1].has_value(),
              scaled_optical_depth[i],
              slant_optical_depths[i + 1] - slant_optical_depths[i],
              LARGE_COSINE_SENTINEL);

          const double BEAM_AT_TOP = std::exp(-slant_optical_depths[i]);
          const double BEAM_AT_BOTTOM = std::exp(-slant_optical_depths[i + 1]);

          const auto COEFF = ComputeLayerCoefficients(
              scaled_asymmetry[i],
              scaled_single_scattering_albedo[i],
              scaled_optical_depth[i],
              COSINE_SOLAR_ZENITH,
              effective_cosine[i + 1],
              BEAM_AT_TOP,
              BEAM_AT_BOTTOM,
              SINGULARITY_GUARD,
              SOURCE_DIVISOR_FLOOR);

          lambda[i] = COEFF.lambda_;
          gamma_ratio[i] = COEFF.gamma_ratio_;
          diffuse_cosine[i] = COEFF.diffuse_cosine_;
          e1[i] = COEFF.e1_;
          e2[i] = COEFF.e2_;
          e3[i] = COEFF.e3_;
          e4[i] = COEFF.e4_;
          c_plus_top[i] = COEFF.c_plus_top_;
          c_minus_top[i] = COEFF.c_minus_top_;
          c_plus_bottom[i] = COEFF.c_plus_bottom_;
          c_minus_bottom[i] = COEFF.c_minus_bottom_;
        }

        // ── Assemble tridiagonal system (Toon et al. eqs. 39–43) ──────────
        const double SURFACE_SOLAR_SOURCE =
            SURFACE_REFLECTIVITY * COSINE_SOLAR_ZENITH * std::exp(-slant_optical_depths[N_LAYERS]);

        TridiagonalMatrix<double> tridiagonal(SYSTEM_SIZE);
        Array1D<double> rhs(SYSTEM_SIZE);

        // Top boundary: zero diffuse incidence (Toon et al. eq. 39)
        tridiagonal.main_diagonal_[0] = e1[0];
        tridiagonal.upper_diagonal_[0] = -e2[0];
        rhs[0] = -c_minus_top[0];

        // Interior rows: flux continuity at layer interfaces
        for (std::size_t k = 0; k < N_LAYERS - 1; ++k)
        {
          const std::size_t EVEN_ROW = (2 * k) + 1;
          tridiagonal.lower_diagonal_[EVEN_ROW - 1] = (e2[k + 1] * e1[k]) - (e3[k] * e4[k + 1]);
          tridiagonal.main_diagonal_[EVEN_ROW] = (e2[k] * e2[k + 1]) - (e4[k] * e4[k + 1]);
          tridiagonal.upper_diagonal_[EVEN_ROW] = (e1[k + 1] * e4[k + 1]) - (e2[k + 1] * e3[k + 1]);
          rhs[EVEN_ROW] =
              ((c_plus_top[k + 1] - c_plus_bottom[k]) * e2[k + 1]) - ((c_minus_top[k + 1] - c_minus_bottom[k]) * e4[k + 1]);

          const std::size_t ODD_ROW = (2 * k) + 2;
          tridiagonal.lower_diagonal_[ODD_ROW - 1] = (e2[k] * e3[k]) - (e4[k] * e1[k]);
          tridiagonal.main_diagonal_[ODD_ROW] = (e1[k] * e1[k + 1]) - (e3[k] * e3[k + 1]);
          tridiagonal.upper_diagonal_[ODD_ROW] = (e3[k] * e4[k + 1]) - (e1[k] * e2[k + 1]);
          rhs[ODD_ROW] =
              (e3[k] * (c_plus_top[k + 1] - c_plus_bottom[k])) + (e1[k] * (c_minus_bottom[k] - c_minus_top[k + 1]));
        }

        // Bottom boundary: surface albedo condition (Toon et al. eq. 43)
        tridiagonal.lower_diagonal_[SYSTEM_SIZE - 2] = e1[N_LAYERS - 1] - (SURFACE_REFLECTIVITY * e3[N_LAYERS - 1]);
        tridiagonal.main_diagonal_[SYSTEM_SIZE - 1] = e2[N_LAYERS - 1] - (SURFACE_REFLECTIVITY * e4[N_LAYERS - 1]);
        rhs[SYSTEM_SIZE - 1] =
            SURFACE_SOLAR_SOURCE - c_plus_bottom[N_LAYERS - 1] + (SURFACE_REFLECTIVITY * c_minus_bottom[N_LAYERS - 1]);

        tuvx::Solve(tridiagonal, rhs);
        const Array1D<double>& solution = rhs;

        // ── Back-substitute: compute fluxes (internal top-to-bottom order) ─
        std::vector<double> direct_irradiance(N_LEVELS);
        std::vector<double> upwelling_irradiance(N_LEVELS);
        std::vector<double> downwelling_irradiance(N_LEVELS);
        std::vector<double> direct_actinic_flux(N_LEVELS);
        std::vector<double> upwelling_actinic_flux(N_LEVELS);
        std::vector<double> downwelling_actinic_flux(N_LEVELS);

        // Level 0 = TOA
        direct_actinic_flux[0] = std::exp(-slant_optical_depths[0]);
        direct_irradiance[0] = COSINE_SOLAR_ZENITH * direct_actinic_flux[0];
        downwelling_irradiance[0] = 0.0;
        upwelling_irradiance[0] = (solution[0] * e3[0]) - (solution[1] * e4[0]) + c_plus_top[0];
        downwelling_actinic_flux[0] = downwelling_irradiance[0] / diffuse_cosine[0];
        upwelling_actinic_flux[0] = upwelling_irradiance[0] / diffuse_cosine[0];

        for (std::size_t lev = 1; lev <= N_LAYERS; ++lev)
        {
          const std::size_t LAYER = lev - 1;
          const std::size_t ROW_IDX = 2 * LAYER;
          direct_actinic_flux[lev] = std::exp(-slant_optical_depths[lev]);
          direct_irradiance[lev] = COSINE_SOLAR_ZENITH * direct_actinic_flux[lev];
          downwelling_irradiance[lev] =
              (solution[ROW_IDX] * e3[LAYER]) + (solution[ROW_IDX + 1] * e4[LAYER]) + c_minus_bottom[LAYER];
          upwelling_irradiance[lev] =
              (solution[ROW_IDX] * e1[LAYER]) + (solution[ROW_IDX + 1] * e2[LAYER]) + c_plus_bottom[LAYER];
          downwelling_actinic_flux[lev] = downwelling_irradiance[lev] / diffuse_cosine[LAYER];
          upwelling_actinic_flux[lev] = upwelling_irradiance[lev] / diffuse_cosine[LAYER];
        }

        // ── Write results reversed to bottom-to-top convention ─────────────
        for (std::size_t lev = 0; lev <= N_LAYERS; ++lev)
        {
          const std::size_t OUT_LEVEL = N_LAYERS - lev;
          radiation_field.spectral_irradiance_.direct_(wl, OUT_LEVEL, col) = direct_irradiance[lev];
          radiation_field.spectral_irradiance_.upwelling_(wl, OUT_LEVEL, col) = upwelling_irradiance[lev];
          radiation_field.spectral_irradiance_.downwelling_(wl, OUT_LEVEL, col) = downwelling_irradiance[lev];
          radiation_field.actinic_flux_.direct_(wl, OUT_LEVEL, col) = direct_actinic_flux[lev];
          radiation_field.actinic_flux_.upwelling_(wl, OUT_LEVEL, col) = upwelling_actinic_flux[lev];
          radiation_field.actinic_flux_.downwelling_(wl, OUT_LEVEL, col) = downwelling_actinic_flux[lev];
        }
      }  // wavelength loop
    }  // column loop
  }

}  // namespace tuvx
