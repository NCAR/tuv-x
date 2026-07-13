// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
//
// Spherical geometry for slant-path optical depth calculations.
#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numbers>
#include <optional>
#include <vector>

namespace tuvx
{

  /// @brief Slant-path geometry calculator for a spherical atmosphere.
  ///
  /// Computes the ratio of slant-path length to vertical layer depth for the
  /// direct solar beam through each atmospheric layer, accounting for Earth's
  /// curvature.
  ///
  /// \rst
  /// Algorithm from :cite:`DahlbackStamnes1991`.
  /// \endrst
  ///
  /// Index conventions:
  ///   - Level index 0 = top of atmosphere (TOA); level n_layers = ground.
  ///   - Layer index 0 = topmost layer (between levels 0 and 1).
  ///   - Altitude edges are ordered bottom-to-top (index 0 = ground).
  ///
  /// All inputs and outputs are in SI units (meters, radians).
  struct SphericalGeometry
  {
    /// Number of layers crossed by the direct beam at each level (0 = TOA).
    /// @c std::nullopt means the level is below the tangent height (no direct beam).
    std::vector<std::optional<std::size_t>> nid_;

    /// Slant-path / vertical-depth ratio for each (level, layer) pair.
    /// @c dsdh_[level][layer] is the ratio for the direct beam reaching @c level
    /// through layer @c layer, with layer indices ordered top-to-bottom (0-based).
    std::vector<std::vector<double>> dsdh_;

    /// @brief Compute slant-path geometry for a given solar zenith angle.
    /// @param solar_zenith_angle Solar zenith angle [radians].
    /// @param altitude_edges Altitude at each grid edge [m], ordered
    ///        bottom-to-top (index 0 = ground, index n_layers = TOA).
    void SetParameters(double solar_zenith_angle, const std::vector<double>& altitude_edges);
  };

  /// @brief Total slant-path optical depth from TOA to a given level.
  ///
  /// @param level Level index (0 = TOA). Layers 0..level-1 lie above this level.
  /// @param n_layers_crossed Layers crossed by the beam (SphericalGeometry::nid_
  ///        at this level). @c std::nullopt → level below tangent height; returns infinity.
  /// @param slant_path Per-layer slant/vertical ratios at this level
  ///        (SphericalGeometry::dsdh_[level]), ordered top-to-bottom, 0-based.
  /// @param optical_depth Scaled layer optical depths, ordered top-to-bottom.
  /// @return Total slant-path optical depth.
  [[nodiscard]] double SlantOpticalDepth(  // NOLINT(readability-identifier-naming)
      std::size_t level,
      std::optional<std::size_t> n_layers_crossed,
      const std::vector<double>& slant_path,
      const std::vector<double>& optical_depth);

  // ── Inline implementations ────────────────────────────────────────────────

  inline void SphericalGeometry::SetParameters(double solar_zenith_angle, const std::vector<double>& altitude_edges)
  {
    static constexpr double EARTH_RADIUS = 6.371e6;
    static constexpr double HALF_PI = std::numbers::pi / 2.0;

    const std::size_t N_LAYERS = altitude_edges.size() - 1;
    const double SURFACE_ELEVATION = altitude_edges[0];
    const double EARTH_RADIUS_AT_SURFACE = EARTH_RADIUS + SURFACE_ELEVATION;
    const bool ABOVE_HORIZON = (solar_zenith_angle <= HALF_PI);
    const double SIN_SZA = std::sin(solar_zenith_angle);

    // Reorder altitudes top-to-bottom: altitude_from_toa[0] = TOA, [n_layers] = ground
    std::vector<double> altitude_from_toa(N_LAYERS + 1);
    for (std::size_t idx = 0; idx <= N_LAYERS; ++idx)
    {
      altitude_from_toa[idx] = altitude_edges[N_LAYERS - idx] - SURFACE_ELEVATION;
    }

    nid_.assign(N_LAYERS + 1, std::nullopt);
    dsdh_.assign(N_LAYERS + 1, std::vector<double>(N_LAYERS, 0.0));

    // Returns the number of layers the direct beam crosses to reach 'level',
    // or nullopt if the level is below the tangent height. Captures shared geometry.
    auto find_layers_crossed = [&](std::size_t level) -> std::optional<std::size_t>
    {
      const double IMPACT_PARAMETER = (EARTH_RADIUS_AT_SURFACE + altitude_from_toa[level]) * SIN_SZA;
      if (!ABOVE_HORIZON && IMPACT_PARAMETER < EARTH_RADIUS_AT_SURFACE)
      {
        return std::nullopt;
      }
      if (ABOVE_HORIZON)
      {
        return level;
      }
      for (std::size_t j = 1; j <= N_LAYERS; ++j)
      {
        if (IMPACT_PARAMETER < altitude_from_toa[j - 1] + EARTH_RADIUS_AT_SURFACE &&
            IMPACT_PARAMETER >= altitude_from_toa[j] + EARTH_RADIUS_AT_SURFACE)
        {
          return j;
        }
      }
      return std::nullopt;
    };

    // Fills dsdh_[level] with slant/vertical-depth ratios. Reads nid_[level].
    auto compute_slant_ratios = [&](std::size_t level)
    {
      const double IMPACT_PARAMETER = (EARTH_RADIUS_AT_SURFACE + altitude_from_toa[level]) * SIN_SZA;
      const std::size_t LAYERS_CROSSED = nid_[level].value();
      for (std::size_t j = 1; j <= LAYERS_CROSSED; ++j)
      {
        const double RADIUS_UPPER = EARTH_RADIUS_AT_SURFACE + altitude_from_toa[j - 1];
        const double RADIUS_LOWER = EARTH_RADIUS_AT_SURFACE + altitude_from_toa[j];
        const double LAYER_THICKNESS = altitude_from_toa[j - 1] - altitude_from_toa[j];
        const double UPPER_TERM = std::max(0.0, (RADIUS_UPPER * RADIUS_UPPER) - (IMPACT_PARAMETER * IMPACT_PARAMETER));
        const double LOWER_TERM = std::max(0.0, (RADIUS_LOWER * RADIUS_LOWER) - (IMPACT_PARAMETER * IMPACT_PARAMETER));
        // The tangent-ray layer contributes negatively when SZA > 90°: the beam
        // enters and exits on the same side of the layer.
        const double PATH_SIGN = (j == LAYERS_CROSSED && LAYERS_CROSSED == level && !ABOVE_HORIZON) ? -1.0 : 1.0;
        const double SLANT_PATH_LENGTH = PATH_SIGN * (std::sqrt(UPPER_TERM) - std::sqrt(LOWER_TERM));
        if (LAYER_THICKNESS > 0.0)
        {
          dsdh_[level][j - 1] = SLANT_PATH_LENGTH / LAYER_THICKNESS;
        }
      }
    };

    for (std::size_t i = 0; i <= N_LAYERS; ++i)
    {
      nid_[i] = find_layers_crossed(i);
      if (nid_[i].has_value())
      {
        compute_slant_ratios(i);
      }
    }
  }

  inline double SlantOpticalDepth(  // NOLINT(readability-identifier-naming)
      std::size_t level,
      std::optional<std::size_t> n_layers_crossed,
      const std::vector<double>& slant_path,
      const std::vector<double>& optical_depth)
  {
    if (!n_layers_crossed.has_value())
    {
      return std::numeric_limits<double>::infinity();
    }
    double result = 0.0;
    const std::size_t DIRECT_PATH_LAYERS = std::min(n_layers_crossed.value(), level);
    for (std::size_t j = 0; j < DIRECT_PATH_LAYERS; ++j)
    {
      result += (optical_depth[j] * slant_path[j]);
    }
    for (std::size_t j = DIRECT_PATH_LAYERS; j < n_layers_crossed.value(); ++j)
    {
      result += 2.0 * (optical_depth[j] * slant_path[j]);
    }
    return result;
  }

}  // namespace tuvx
