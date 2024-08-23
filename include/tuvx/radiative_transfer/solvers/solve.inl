#include <tuvx/radiative_transfer/solvers/delta_eddington.hpp>

// Copyright (C) 2023-2024 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

namespace tuvx
{
  template<
      typename T,
      typename ArrayPolicy,
      typename GridPolicy,
      typename ProfilePolicy,
      typename RadiatorStatePolicy,
      typename RadiationFieldPolicy>
  inline void Solve(
      const std::vector<T>& solar_zenith_angles,
      const std::map<std::string, GridPolicy>& grids,
      const std::map<std::string, ProfilePolicy>& profiles,
      const std::function<void(const RadiatorStatePolicy&, const ArrayPolicy&, const std::vector<T>)> ApproximationFunction,
      const RadiatorStatePolicy& accumulated_radiator_states,
      RadiationFieldPolicy& radiation_field){
    // Solve the radiative transfer equation.
    //
    // Things that will change from the original solver:
    // 1. All variables will be in SI units. Some of the original solver's
    //    variables were in non-SI units.
    // 2. We will be solving for collections of columns. The original solver
    //    was for a single column.
    // 3. The variable naming and source-code documentation will be improved.

    const std::size_t number_of_columns = solar_zenith_angles.size();
    const auto& vertical_grid = grids.at("altitude [m]");
    const auto& wavelength_grid = grids.at("wavelength [m]");

    // Check for consistency between the grids and profiles.
    assert(vertical_grid.NumberOfColumns() == number_of_columns);
    assert(wavelength_grid.NumberOfColumns() == 1);

    // dictionaryu to hold internal solver variables
    std::map<std::string, ArrayPolicy> solver_variables;

    // dictionary to hold source functions (C1, C2 functions from the paper)
    std::map<std::string, Array3D<RadiatorStatePolicy>> source_functions;

    // Initlialize variables (delta scaling, source functions etc)
    InitializeVariables<T, GridPolicy, ProfilePolicy, RadiatorStatePolicy>(
        solar_zenith_angles, grids, profiles, accumulated_radiator_states, solver_variables, source_functions);
  }
}  // namespace tuvx
