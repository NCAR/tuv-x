// Copyright (C) 2023-2024 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

#include <tuvx/radiative_transfer/solvers/delta_eddington.hpp>

namespace tuvx {

  template <typename T>
  RadiationField<T> DeltaEddington::Solve(const std::vector<T>& solar_zenith_angles,
                             const std::map<std::string, Grid<T>>& grids,
                             const std::map<std::string, Profile<T>>& profiles,
                             const RadiatorState<T>& accumulated_radiator_state) const {
    // Solve the radiative transfer equation.
    // 
    // [DEV NOTES] This is a placeholder for the actual implementation.
    // The spherical geometry argument of the original solver was left out
    // until we determine whether it needs to be an object or just a set of functions.
    //
    // Things that will change from the original solver:
    // 1. All variables will be in SI units. Some of the original solver's
    //    variables were in non-SI units.
    // 2. We will be solving for collections of columns. The original solver
    //    was for a single column.
    const size_t number_of_columns = solar_zenith_angles.size();
    const auto& vertical_grid = grids.at("altitude [m]");
    const auto& wavelength_grid = grids.at("wavelength [m]");

    // Check for consistency between the grids and profiles.
    assert(vertical_grid.number_of_columns() == number_of_columns);
    assert(wavelength_grid.number_of_columns() == 1);

    RadiationField<T> radiation_field(number_of_columns, vertical_grid, wavelength_grid);

    return radiation_field;
  }

} // namespace tuvx