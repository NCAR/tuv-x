// Copyright (C) 2023-2024 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
//
/// @file Radiation field vertically and wavelength resolved.
#pragma once

#include <vector>

#include <tuvx/util/array3d.hpp>
#include <tuvx/grid.hpp>

namespace tuvx {

  /// Radiation field components vertically and wavelength resolved.
  ///
  /// Data structures are designed to hold the radiation field components
  /// for a collection of columns in a 3D grid.
  /// Arays are indexed by column, vertical edge, and wavelength.
  template <typename T=double>
  struct RadiationFieldComponents {
    /// Direct component of the radiation field.
    Array3D<T> direct_;
    /// Upwelling component of the radiation field.
    Array3D<T> upwelling_;
    /// Downwelling component of the radiation field.
    Array3D<T> downwelling_;

    /// Constructor
    RadiationFieldComponents(size_t number_of_columns, Grid<T>& vertical_grid,
                             Grid<T>& wavelength_grid)
      : direct_(number_of_columns, vertical_grid.number_of_edges(),
                wavelength_grid.number_of_cells()),
        upwelling_(number_of_columns, vertical_grid.number_of_edges(),
                   wavelength_grid.number_of_cells()),
        downwelling_(number_of_columns, vertical_grid.number_of_edges(),
                     wavelength_grid.number_of_cells()) {}
  };

  /// Radiation field vertically and wavelength resolved.
  template <typename T=double>
  struct RadiationField {
    /// Total spectral irradiance.
    RadiationFieldComponents<T> spectral_irradiance_;
    /// Total actinic flux.
    RadiationFieldComponents<T> actinic_flux_;

    /// Constructor
    RadiationField(size_t number_of_columns, Grid<T>& vertical_grid,
                   Grid<T>& wavelength_grid)
      : spectral_irradiance_(number_of_columns, vertical_grid, wavelength_grid),
        actinic_flux_(number_of_columns, vertical_grid, wavelength_grid) {}
  };

} // namespace tuvx