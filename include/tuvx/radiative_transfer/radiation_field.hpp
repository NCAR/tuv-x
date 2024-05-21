// Copyright (C) 2023-2024 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
//
/// @file Radiation field vertically and wavelength resolved.
#pragma once

#include <vector>

#include <tuvx/util/array3d.hpp>
#include <tuvx/grid.hpp>

namespace tuvx {

  /// @brief Radiation field components vertically and wavelength resolved.
  ///
  /// Data structures are designed to hold the radiation field components
  /// for a collection of columns in a 3D grid.
  /// Arays are indexed by column, vertical edge, and wavelength.
  template <typename ArrayPolicy = Array3D<double>>
  struct RadiationFieldComponents {
    /// Direct component of the radiation field.
    ArrayPolicy direct_;
    /// Upwelling component of the radiation field.
    ArrayPolicy upwelling_;
    /// Downwelling component of the radiation field.
    ArrayPolicy downwelling_;

    /// @brief Constructor for radiation field components.
    /// @param number_of_columns Number of columns.
    /// @param vertical_grid Vertical grid.
    /// @param wavelength_grid Wavelength grid.
    template<typename GridPolicy>
    RadiationFieldComponents(size_t number_of_columns, GridPolicy& vertical_grid,
                             GridPolicy& wavelength_grid)
      : direct_(number_of_columns, vertical_grid.number_of_edges(),
                wavelength_grid.number_of_cells()),
        upwelling_(number_of_columns, vertical_grid.number_of_edges(),
                   wavelength_grid.number_of_cells()),
        downwelling_(number_of_columns, vertical_grid.number_of_edges(),
                     wavelength_grid.number_of_cells()) {}
  };

  /// @brief Radiation field vertically and wavelength resolved.
  template <typename ComponentPolicy = RadiationFieldComponents<Array3D<double>>>
  struct RadiationField {
    /// Total spectral irradiance.
    ComponentPolicy spectral_irradiance_;
    /// Total actinic flux.
    ComponentPolicy actinic_flux_;

    /// @brief Constructor for a radiation field.
    /// @param number_of_columns Number of columns.
    /// @param vertical_grid Vertical grid.
    /// @param wavelength_grid Wavelength grid.
    template<typename GridPolicy>
    RadiationField(size_t number_of_columns, GridPolicy& vertical_grid,
                   GridPolicy& wavelength_grid)
      : spectral_irradiance_(number_of_columns, vertical_grid, wavelength_grid),
        actinic_flux_(number_of_columns, vertical_grid, wavelength_grid) {}
  };

} // namespace tuvx