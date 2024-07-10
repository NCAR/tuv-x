// Copyright (C) 2023-2024 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#include "tuvx/linear_algebra/linear_algebra.hpp"

#include <tuvx/radiative_transfer/solvers/delta_eddington.hpp>
namespace tuvx
{

  /// @brief Initilize the variables required for the delta eddington approximation
  /// @param solar zenith angles Solar zenith angles for each column [radians]
  /// @param grids Grids available for the radiative transfer calculation
  /// @param profiles Profiles available for the radiative transfer calculation
  /// @param radiation_field The calculated radiation field
  ///
  /// Delta scaling of the inputs, compuation of reflectivity, diffuse flux
  template<
      typename T,
      typename GridPolicy,
      typename ProfilePolicy,
      typename RadiatorStatePolicy,
      typename RadiationFieldPolicy,
      typename ArrayPolicy>
  inline void InitializeVariables(
      const std::vector<T>& solar_zenith_angles,
      const std::map<std::string, GridPolicy>& grids,
      const std::map<std::string, ProfilePolicy>& profiles,
      const RadiatorStatePolicy& accumulated_radiator_states)
  {
    // determine number of layers
    const std::size_t number_of_columns = solar_zenith_angles.size();
    const auto& vertical_grid = grids.at("altitude [m]");
    const auto& wavelength_grid = grids.at("wavelength [m]");

    // Check for consistency between the grids and profiles.
    assert(vertical_grid.NumberOfColumns() == number_of_columns);
    assert(wavelength_grid.NumberOfColumns() == 1);

    // Radiator state variables
    ArrayPolicy& optical_depth = accumulated_radiator_states.optical_depth_;
    ArrayPolicy& single_scattering_albedo = accumulated_radiator_states.single_scattering_albedo_;
    ArrayPolicy& assymetry_parameter = accumulated_radiator_states.assymetry_parameter_;

    // Delta scaling
    T f;
    for (std::size_t i = 0; i < number_of_columns; i++)
    {
      f = assymetry_parameter[i] * assymetry_parameter[i];
      assymetry_parameter[i] = (assymetry_parameter - f) / (1 - f);
      single_scattering_albedo[i] = (1 - f) * single_scattering_albedo[i] / (1 - single_scattering_albedo[i] * f);
      optical_depth[i] = (1 - single_scattering_albedo[i] * f) * optical_depth[i];
    }

    // calculate slant optical depth
    for (auto& tau : optical_depth)
    {
      // TODO slant optical depth computation
      tau = tau;
    }

    // delta eddington parameters (little gammas 1, 2, 3 and mu from the paper)
    ArrayPolicy delta_eddington_parameters(number_of_columns, 5);
    {
      for (std::size_t i = 0; i < number_of_columns; i++)
      {
        // delta eddington parameters
        auto& gamma1 = delta_eddington_parameters[i][0];
        auto& gamma2 = delta_eddington_parameters[i][1];
        auto& gamma3 = delta_eddington_parameters[i][2];
        auto& gamma4 = delta_eddington_parameters[i][3];
        auto& mu = delta_eddington_parameters[i][4];

        // radiator state variables
        auto mu_0 = std::acos(solar_zenith_angles[i]);
        auto& omega = optical_depth[i];
        auto& g = single_scattering_albedo[i];
        auto& tau = optical_depth[i];

        // compute delta eddington parameters
        gamma1 = 7 - omega * (4 + 3 * g);
        gamma2 = -(1 - omega * (4 - 3 * g)) / 4;
        gamma3 = (2 - 3 * g * mu_0) / 4;
        mu = (T)0.5;
      }
    }
  }

  template<
      typename T,
      typename GridPolicy,
      typename ProfilePolicy,
      typename RadiatorStatePolicy,
      typename RadiationFieldPolicy>
  inline void DeltaEddington::Solve(
      const std::vector<T>& solar_zenith_angles,
      const std::map<std::string, GridPolicy>& grids,
      const std::map<std::string, ProfilePolicy>& profiles,
      const RadiatorStatePolicy& accumulated_radiator_state,
      RadiationFieldPolicy& radiation_field) const
  {
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
    const std::size_t number_of_columns = solar_zenith_angles.size();
    const auto& vertical_grid = grids.at("altitude [m]");
    const auto& wavelength_grid = grids.at("wavelength [m]");

    // Check for consistency between the grids and profiles.
    assert(vertical_grid.NumberOfColumns() == number_of_columns);
    assert(wavelength_grid.NumberOfColumns() == 1);

    // internal solver variables
    Array2D<T> solution_parameters;
    Array2D<T> delta_eddington_parameters;

    // tridiagonal system variables
    TridiagonalMatrix<T> coeffcient_matrix;
    std::vector<T> coeffcient_vector;

    // Initialize variables
    this->InitializeVariables<T, GridPolicy, ProfilePolicy, RadiatorStatePolicy, RadiationFieldPolicy>(
        solar_zenith_angles, grids, profiles, accumulated_radiator_state);

    // Assemble the tridiagonal system
    this->AssembleTridiagonalSystem<T, GridPolicy, ProfilePolicy, RadiatorStatePolicy, RadiationFieldPolicy>(
        solar_zenith_angles, grids, profiles, solution_parameters, coeffcient_matrix, coeffcient_vector);

    // solve tridiagonal system
    Solve<T>(coeffcient_matrix, coeffcient_vector);

    // Update radiation field
    this->ComputeRadiationField<T, GridPolicy, ProfilePolicy, RadiatorStatePolicy>(
        solar_zenith_angles, grids, profiles, solution_parameters, coeffcient_matrix, coeffcient_vector, radiation_field);

    // [DEV NOTES] Temporarily return predictable values for the radiation field.
    // This will be replaced with the actual results once the solver is implemented.
    int offset = 42;
    for (auto& elem : radiation_field.spectral_irradiance_.direct_)
    {
      elem = offset++;
    }
    offset = 93;
    for (auto& elem : radiation_field.spectral_irradiance_.upwelling_)
    {
      elem = offset++;
    }
    offset = 52;
    for (auto& elem : radiation_field.spectral_irradiance_.downwelling_)
    {
      elem = offset++;
    }
    offset = 5;
    for (auto& elem : radiation_field.actinic_flux_.direct_)
    {
      elem = offset++;
    }
    offset = 24;
    for (auto& elem : radiation_field.actinic_flux_.upwelling_)
    {
      elem = offset++;
    }
    offset = 97;
    for (auto& elem : radiation_field.actinic_flux_.downwelling_)
    {
      elem = offset++;
    }
  }

}  // namespace tuvx
