// Copyright (C) 2023-2024 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
//
/// @file Optically active atmospheric component
#pragma once

#include <tuvx/util/array3d.hpp>
#include <tuvx/grid.hpp>

namespace tuvx {

  /// Optical properties of a radiative transfer component.
  ///
  /// The optical properties are represented by 3D arrays indexed by
  /// [column][wavelength][vertical layer].
  template <typename T=double>
  class RadiatorState {
  public:
    /// Layer optical depth (unitless)
    Array3D<T> optical_depth_;
    /// Single scattering albedo (unitless)
    Array3D<T> single_scattering_albedo_;
    /// Asymmetry parameter (unitless)
    Array3D<T> asymmetry_parameter_;

    /// Constructor
    RadiatorState(size_t number_of_columns, Grid<T>& vertical_grid, Grid<T>& wavelength_grid)
      : optical_depth_(number_of_columns, wavelength_grid.number_of_cells(),
                       vertical_grid.number_of_cells()),
        single_scattering_albedo_(number_of_columns, wavelength_grid.number_of_cells(),
                                  vertical_grid.number_of_cells()),
        asymmetry_parameter_(number_of_columns, wavelength_grid.number_of_cells(),
                             vertical_grid.number_of_cells()) {}
    
    /// Accumulate a set of optical properties.
    static RadiatorState<T> Accumulate(const std::vector<RadiatorState<T>>& states) {
      // [DEV NOTES] Placeholder for the RadiatorState::Accumulate method
      return states[0];
    }
  };

  // Placeholder for the Radiator class (if needed)

} // namespace tuvx