// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
//
// Atmospheric state container passed to every TransformFunc.
#pragma once

#include <tuvx/grid.hpp>
#include <tuvx/util/array2d.hpp>
#include <tuvx/util/array3d.hpp>

namespace tuvx
{

  /// @brief Atmospheric state used to calculate transform weights.
  ///
  /// Holds the grids and per-column profiles that transforms need to
  /// calculate their weight arrays [wavelength x height x column].
  /// The template parameter controls the 3D storage policy; 2D profiles
  /// are always backed by Array2D<value_type>.
  ///
  /// @tparam ArrayPolicy  3D array policy (default: Array3D<double>).
  template<typename ArrayPolicy = Array3D<double>>
  struct AtmosphericState
  {
    using value_type = typename ArrayPolicy::value_type;

    /// Wavelength grid (constant across columns), [n_wavelengths x 1].
    Grid<Array2D<value_type>> wavelength_grid_{};

    /// Height grid (may vary by column), [n_layers x n_columns].
    Grid<Array2D<value_type>> height_grid_{};

    /// Temperature profile [n_layers x n_columns] (K).
    Array2D<value_type> temperature_{};

    /// Pressure profile [n_layers x n_columns] (Pa).
    Array2D<value_type> pressure_{};

    /// Air number density [n_layers x n_columns] (mol m-3).
    Array2D<value_type> air_density_{};
  };

}  // namespace tuvx
