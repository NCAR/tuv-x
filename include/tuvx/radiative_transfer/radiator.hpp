// Copyright (C) 2023-2024 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
//
/// @file Optically active atmospheric component
#pragma once

#include <tuvx/util/array3d.hpp>
#include <tuvx/grid.hpp>

namespace tuvx {

  /// @brief Optical properties of a radiative transfer component.
  ///
  /// The optical properties are represented by 3D arrays indexed by
  /// [column][wavelength][vertical layer].
  template <typename ArrayPolicy = Array3D<double>>
  class RadiatorState {
  public:
    /// Layer optical depth (unitless)
    ArrayPolicy optical_depth_;
    /// Single scattering albedo (unitless)
    ArrayPolicy single_scattering_albedo_;
    /// Asymmetry parameter (unitless)
    ArrayPolicy asymmetry_parameter_;

    /// @brief Constructor of a radiator state
    /// @param number_of_columns Number of columns.
    /// @param vertical_grid Vertical grid.
    /// @param wavelength_grid Wavelength grid.
    template<typename GridPolicy>
    RadiatorState(size_t number_of_columns, GridPolicy& vertical_grid, GridPolicy& wavelength_grid)
      : optical_depth_(number_of_columns, wavelength_grid.number_of_cells(),
                       vertical_grid.number_of_cells()),
        single_scattering_albedo_(number_of_columns, wavelength_grid.number_of_cells(),
                                  vertical_grid.number_of_cells()),
        asymmetry_parameter_(number_of_columns, wavelength_grid.number_of_cells(),
                             vertical_grid.number_of_cells()) {}
    
    /// @brief Accumulate a set of optical properties.
    /// @param states Vector of optical properties for individual optically active components.
    /// @return Accumulated optical properties.
    template<typename RadiatorStatePolicy>
    static RadiatorStatePolicy Accumulate(const std::vector<RadiatorStatePolicy>& states) {
      // [DEV NOTES] Placeholder for the RadiatorState::Accumulate method
      return states[0];
    }
  };

  // Placeholder for the Radiator class (if needed)

} // namespace tuvx