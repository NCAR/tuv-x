// Copyright (C) 2023-2025 University Corporation for Atmospheric Research-National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
//
// Delta-Eddington solver for radiative transfer
#pragma once

#include <tuvx/grid.hpp>
#include <tuvx/profile.hpp>
#include <tuvx/radiative_transfer/radiation_field.hpp>
#include <tuvx/radiative_transfer/radiator.hpp>

#include <cassert>
#include <map>
#include <vector>

namespace tuvx
{

  /// @brief Radiative flux calculator that applies the delta-Eddington Approximation.
  ///
  /// [DEV NOTES] We can determine whether this should be a class or a set of functions
  class DeltaEddington
  {
   public:
    /// Construct a Delta-Eddington solver.
    DeltaEddington() = default;

    /// @brief Solve the radiative transfer equation for a collection of columns
    /// @param solar_zenith_angles Solar zenith angles for each column [radians].
    /// @param grids Grids available for the radiative transfer calculation.
    /// @param profiles Profiles available for the radiative transfer calculation.
    /// @param radiation_field The calculated radiation field.
    ///
    /// Solves two-stream equations for multiple layers. These routines are based
    /// on equations from: Toon et al., J.Geophys.Res., v94 (D13), Nov 20, 1989.
    /// DOI: https://doi.org/10.1029/JD094iD13p16287
    /// It contains 9 two-stream methods to choose from. A pseudo-spherical
    /// correction has also been added.
    ///
    /// The original delta-Eddington paper is:
    /// Joseph and Wiscombe, J. Atmos. Sci., 33, 2453-2459, 1976
    /// DOI: https://doi.org/10.1175/1520-0469(1976)033%3C2452:TDEAFR%3E2.0.CO;2
    template<
        typename T,
        typename GridPolicy,
        typename ProfilePolicy,
        typename RadiatorStatePolicy,
        typename RadiationFieldPolicy>
    void Solve(
        const std::vector<T>& solar_zenith_angles,
        const std::map<std::string, GridPolicy>& grids,
        const std::map<std::string, ProfilePolicy>& profiles,
        const RadiatorStatePolicy& accumulated_radiator_states,
        RadiationFieldPolicy& radiation_field) const;
  };

}  // namespace tuvx

#include "delta_eddington.inl"