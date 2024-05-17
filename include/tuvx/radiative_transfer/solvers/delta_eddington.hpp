// Copyright (C) 2023-2024 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
//
/// @file Delta-Eddington solver for radiative transfer.
#pragma once

#include <vector>
#include <map>

#include <tuvx/grid.hpp>
#include <tuvx/profile.hpp>
#include <tuvx/radiative_transfer/radiation_field.hpp>
#include <tuvx/radiative_transfer/radiator.hpp>

namespace tuvx {

  /// Radiative flux calculator that applies the delta-Eddington Approximation.
  ///
  /// [DEV NOTES] We can determine whether this should be a class or a set of functions
  class DeltaEddington {
  public:
    /// Construct a Delta-Eddington solver.
    DeltaEddington() = default;

    /// @brief Solve the radiative transfer equation for a collection of columns
    /// @param solar_zenith_angles Solar zenith angles for each column [radians].
    /// @param grids Grids available for the radiative transfer calculation.
    /// @param profiles Profiles available for the radiative transfer calculation.
    /// @return The calculated radiation field.
    ///
    /// Solves two-stream equations for multiple layers. These routines are based
    /// on equations from: Toon et al., J.Geophys.Res., v94 (D13), Nov 20, 1989.
    /// It contains 9 two-stream methods to choose from. A pseudo-spherical
    /// correction has also been added.
    ///
    /// The original delta-Eddington paper is:
    /// Joseph and Wiscombe, J. Atmos. Sci., 33, 2453-2459, 1976
    template <typename T=double>
    RadiationField<T> Solve(const std::vector<T>& solar_zenith_angles,
               const std::map<std::string, Grid<T>>& grids,
               const std::map<std::string, Profile<T>>& profiles,
               const RadiatorState<T>& accumulated_radiator_states) const;
  };

} // namespace tuvx