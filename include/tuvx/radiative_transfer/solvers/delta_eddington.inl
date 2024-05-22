// Copyright (C) 2023-2024 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

namespace tuvx {

  template <typename T, typename GridPolicy, typename ProfilePolicy, typename RadiatorStatePolicy, typename RadiationFieldPolicy>
  inline void DeltaEddington::Solve(const std::vector<T>& solar_zenith_angles,
                             const std::map<std::string, GridPolicy>& grids,
                             const std::map<std::string, ProfilePolicy>& profiles,
                             const RadiatorStatePolicy& accumulated_radiator_state,
                             RadiationFieldPolicy& radiation_field) const {
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
    // 3. The variable naming and source-code documentation will be improved.
    const size_t number_of_columns = solar_zenith_angles.size();
    const auto& vertical_grid = grids.at("altitude [m]");
    const auto& wavelength_grid = grids.at("wavelength [m]");

    // Check for consistency between the grids and profiles.
    assert(vertical_grid.NumberOfColumns() == number_of_columns);
    assert(wavelength_grid.NumberOfColumns() == 1);

  }

} // namespace tuvx