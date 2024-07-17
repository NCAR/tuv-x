

// Copyright (C) 2023-2024 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#include "tuvx/linear_algebra/linear_algebra.hpp"
#include "tuvx/radiative_transfer/radiator.hpp"

#include <tuvx/radiative_transfer/solvers/delta_eddington.hpp>

#include <cmath>

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
      const RadiatorStatePolicy& accumulated_radiator_state,
      RadiationFieldPolicy& radiation_field)
  {
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

    // tridiagonal system variables
    TridiagonalMatrix<T> coeffcient_matrix(number_of_columns, 0);
    std::vector<T> coeffcient_vector(number_of_columns, 0);

    // internal solver variables
    std::map<std::string, std::vector<T>> solution_parameters;
    std::map<std::string, std::vector<T>> simulation_parameters;
    std::map<std::string, std::function<T(T)>> source_functions;
    for (std::size_t i = 0; i < wavelength_grid.NumberOfColumns(); i++)
    {
      solution_parameters.at("Single Scattering Albedo") = accumulated_radiator_state.single_scattering_albedo_[0][i];
      solution_parameters.at("Optical Depth") = accumulated_radiator_state.optical_depth_[0][i];
      solution_parameters.at("Assymetry Parameter") = accumulated_radiator_state.assymetry_parameter_[0][i];
    }

    tuvx::InitializeVariables<T, GridPolicy, ProfilePolicy, RadiatorStatePolicy, RadiationFieldPolicy>(
        solar_zenith_angles, grids, profiles, accumulated_radiator_state);

    // delta eddington approximation or quadrature
    ApproximationFunction(accumulated_radiator_state, solar_zenith_angles, simulation_parameters);

    tuvx::AssembleTridiagonalMatrix<T>(
        solar_zenith_angles, grids, profiles, solution_parameters, coeffcient_matrix, coeffcient_vector);

    tuvx::Solve<T>(coeffcient_matrix, coeffcient_vector);

    tuvx::ComputeRadiationField<T, GridPolicy, ProfilePolicy, RadiatorStatePolicy>(
        solar_zenith_angles, grids, profiles, solution_parameters, coeffcient_matrix, coeffcient_vector, radiation_field);
  }
}  // namespace tuvx
