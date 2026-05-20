// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
//
// Spherical geometry for slant-path optical depth calculations.
// Ported from src/spherical_geometry.F90 and src/radiative_transfer/solver.F90.
#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <vector>

namespace tuvx
{

  /// @brief Slant-path geometry calculator for a spherical atmosphere.
  ///
  /// Computes the ratio of slant-path length to vertical layer depth for the
  /// direct solar beam through each atmospheric layer, accounting for Earth's
  /// curvature. Based on Dahlback and Stamnes (1991),
  /// doi:10.1016/0032-0633(91)90061-E.
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
    /// A value < 0 means the level is below the tangent height (no direct beam).
    std::vector<int> nid_;

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
  ///        at this level). Negative → level below tangent height; returns 1e36.
  /// @param slant_path Per-layer slant/vertical ratios at this level
  ///        (SphericalGeometry::dsdh_[level]), ordered top-to-bottom, 0-based.
  /// @param optical_depth Scaled layer optical depths, ordered top-to-bottom.
  /// @return Total slant-path optical depth.
  double SlantOpticalDepth(
      std::size_t level,
      int n_layers_crossed,
      const std::vector<double>& slant_path,
      const std::vector<double>& optical_depth);

  // ── Inline implementations ────────────────────────────────────────────────

  inline void SphericalGeometry::SetParameters(
      double solar_zenith_angle,
      const std::vector<double>& altitude_edges)
  {
    // Earth radius in SI units (Fortran uses 6.371e3 km)
    static constexpr double kEarthRadius = 6.371e6;
    static constexpr double kHalfPi      = M_PI / 2.0;

    const std::size_t n_layers       = altitude_edges.size() - 1;
    const double      surface_elev   = altitude_edges[0];
    const double      re             = kEarthRadius + surface_elev;
    const bool        above_horizon  = (solar_zenith_angle <= kHalfPi);
    const double      sin_sza        = std::sin(solar_zenith_angle);

    // Inverted altitude: zd[0] = TOA, zd[n_layers] = 0 (ground level)
    std::vector<double> zd(n_layers + 1);
    for (std::size_t k = 0; k <= n_layers; ++k)
    {
      zd[k] = altitude_edges[n_layers - k] - surface_elev;
    }

    nid_.assign(n_layers + 1, 0);
    dsdh_.assign(n_layers + 1, std::vector<double>(n_layers, 0.0));

    for (std::size_t i = 0; i <= n_layers; ++i)
    {
      const double rpsinz = (re + zd[i]) * sin_sza;

      int id = 0;
      if (!above_horizon && rpsinz < re)
      {
        id = -1;
      }
      else
      {
        id = static_cast<int>(i);
        if (!above_horizon)
        {
          id = -1;
          for (std::size_t j = 1; j <= n_layers; ++j)
          {
            if (rpsinz < zd[j - 1] + re && rpsinz >= zd[j] + re)
            {
              id = static_cast<int>(j);
            }
          }
        }
        for (int j = 1; j <= id; ++j)
        {
          double       sm   = 1.0;
          const double rj   = re + zd[j - 1];
          const double rjp1 = re + zd[static_cast<std::size_t>(j)];
          const double dhj  = zd[j - 1] - zd[static_cast<std::size_t>(j)];
          const double ga   = std::max(0.0, rj * rj - rpsinz * rpsinz);
          const double gb   = std::max(0.0, rjp1 * rjp1 - rpsinz * rpsinz);

          // For the tangent-ray layer when zenith > 90°, the beam enters and
          // exits on the same side, so the path contribution is subtracted.
          if (j == id && id == static_cast<int>(i) && !above_horizon)
          {
            sm = -1.0;
          }

          const double dsj = sm * (std::sqrt(ga) - std::sqrt(gb));
          if (dhj > 0.0)
          {
            dsdh_[i][static_cast<std::size_t>(j) - 1] = dsj / dhj;
          }
        }
      }
      nid_[i] = id;
    }
  }

  inline double SlantOpticalDepth(
      std::size_t               level,
      int                       n_layers_crossed,
      const std::vector<double>& slant_path,
      const std::vector<double>& optical_depth)
  {
    static constexpr double kInfinity = 1.0e36;
    if (n_layers_crossed < 0)
    {
      return kInfinity;
    }
    double    result    = 0.0;
    const int min_count = std::min(n_layers_crossed, static_cast<int>(level));
    for (int j = 0; j < min_count; ++j)
    {
      result += optical_depth[static_cast<std::size_t>(j)] * slant_path[static_cast<std::size_t>(j)];
    }
    for (int j = min_count; j < n_layers_crossed; ++j)
    {
      result +=
          2.0 * optical_depth[static_cast<std::size_t>(j)] * slant_path[static_cast<std::size_t>(j)];
    }
    return result;
  }

}  // namespace tuvx
