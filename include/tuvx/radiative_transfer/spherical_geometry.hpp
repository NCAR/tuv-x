// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
//
// Spherical geometry for slant-path optical depth calculations.
// Ported from src/spherical_geometry.F90 and src/radiative_transfer/solver.F90.
#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
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
  /// Ported from @c slant_optical_depth() in @c src/radiative_transfer/solver.F90.
  ///
  /// @param level Level index (0 = TOA). Layers 0..level-1 lie above this level.
  /// @param n_layers_crossed Layers crossed by the beam (SphericalGeometry::nid_
  ///        at this level). @c std::nullopt → level below tangent height; returns infinity.
  /// @param slant_path Per-layer slant/vertical ratios at this level
  ///        (SphericalGeometry::dsdh_[level]), ordered top-to-bottom, 0-based.
  /// @param optical_depth Scaled layer optical depths, ordered top-to-bottom.
  /// @return Total slant-path optical depth.
  double SlantOpticalDepth(
      std::size_t level,
      std::optional<std::size_t> n_layers_crossed,
      const std::vector<double>& slant_path,
      const std::vector<double>& optical_depth);

  // ── Inline implementations ────────────────────────────────────────────────

  inline void SphericalGeometry::SetParameters(
      double solar_zenith_angle,
      const std::vector<double>& altitude_edges)
  {
    static constexpr double earth_radius = 6.371e6;
    static constexpr double half_pi      = M_PI / 2.0;

    const std::size_t n_layers               = altitude_edges.size() - 1;
    const double      surface_elevation      = altitude_edges[0];
    const double      earth_radius_at_surface = earth_radius + surface_elevation;
    const bool        above_horizon          = (solar_zenith_angle <= half_pi);
    const double      sin_sza               = std::sin(solar_zenith_angle);

    // Reorder altitudes top-to-bottom: altitude_from_toa[0] = TOA, [n_layers] = ground
    std::vector<double> altitude_from_toa(n_layers + 1);
    for (std::size_t idx = 0; idx <= n_layers; ++idx)
    {
      altitude_from_toa[idx] = altitude_edges[n_layers - idx] - surface_elevation;
    }

    nid_.assign(n_layers + 1, std::nullopt);
    dsdh_.assign(n_layers + 1, std::vector<double>(n_layers, 0.0));

    // Returns the number of layers the direct beam crosses to reach 'level',
    // or nullopt if the level is below the tangent height. Captures shared geometry.
    auto find_layers_crossed = [&](std::size_t level) -> std::optional<std::size_t>
    {
      const double impact_parameter = (earth_radius_at_surface + altitude_from_toa[level]) * sin_sza;
      if (!above_horizon && impact_parameter < earth_radius_at_surface)
      {
        return std::nullopt;
      }
      if (above_horizon)
      {
        return level;
      }
      for (std::size_t j = 1; j <= n_layers; ++j)
      {
        if (impact_parameter <  altitude_from_toa[j - 1] + earth_radius_at_surface &&
            impact_parameter >= altitude_from_toa[j]     + earth_radius_at_surface)
        {
          return j;
        }
      }
      return std::nullopt;
    };

    // Fills dsdh_[level] with slant/vertical-depth ratios. Reads nid_[level].
    auto compute_slant_ratios = [&](std::size_t level)
    {
      const double      impact_parameter = (earth_radius_at_surface + altitude_from_toa[level]) * sin_sza;
      const std::size_t layers_crossed   = nid_[level].value();
      for (std::size_t j = 1; j <= layers_crossed; ++j)
      {
        const double radius_upper    = earth_radius_at_surface + altitude_from_toa[j - 1];
        const double radius_lower    = earth_radius_at_surface + altitude_from_toa[j];
        const double layer_thickness = altitude_from_toa[j - 1] - altitude_from_toa[j];
        const double upper_term      = std::max(0.0, radius_upper * radius_upper - impact_parameter * impact_parameter);
        const double lower_term      = std::max(0.0, radius_lower * radius_lower - impact_parameter * impact_parameter);
        // The tangent-ray layer contributes negatively when SZA > 90°: the beam
        // enters and exits on the same side of the layer.
        const double path_sign =
            (j == layers_crossed && layers_crossed == level && !above_horizon)
            ? -1.0 : 1.0;
        const double slant_path_length = path_sign * (std::sqrt(upper_term) - std::sqrt(lower_term));
        if (layer_thickness > 0.0)
        {
          dsdh_[level][j - 1] = slant_path_length / layer_thickness;
        }
      }
    };

    for (std::size_t i = 0; i <= n_layers; ++i)
    {
      nid_[i] = find_layers_crossed(i);
      if (nid_[i].has_value())
      {
        compute_slant_ratios(i);
      }
    }
  }

  inline double SlantOpticalDepth(
      std::size_t                level,
      std::optional<std::size_t> n_layers_crossed,
      const std::vector<double>& slant_path,
      const std::vector<double>& optical_depth)
  {
    if (!n_layers_crossed.has_value())
    {
      return std::numeric_limits<double>::infinity();
    }
    double            result             = 0.0;
    const std::size_t direct_path_layers = std::min(n_layers_crossed.value(), level);
    for (std::size_t j = 0; j < direct_path_layers; ++j)
    {
      result += optical_depth[j] * slant_path[j];
    }
    for (std::size_t j = direct_path_layers; j < n_layers_crossed.value(); ++j)
    {
      result += 2.0 * optical_depth[j] * slant_path[j];
    }
    return result;
  }

}  // namespace tuvx
