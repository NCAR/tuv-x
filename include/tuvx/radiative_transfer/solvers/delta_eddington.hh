// Copyright (C) 2023-2024 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#include "tuvx/linear_algebra/linear_algebra.hpp"
#include "tuvx/radiative_transfer/radiator.hpp"

#include <tuvx/radiative_transfer/solvers/delta_eddington.hpp>

template<typename T, typename RadiatorStatePolicy, typename ArrayPolicy>
void DeltaEddingtonApproximation(
    const RadiatorStatePolicy& accumulated_radiator_states,
    const ArrayPolicy& parameters,
    const std::vector<T> solar_zenith_angles)
{
  const std::size_t number_of_columns = solar_zenith_angles.size();
  ArrayPolicy& optical_depth = accumulated_radiator_states.optical_depth_;
  ArrayPolicy& single_scattering_albedo = accumulated_radiator_states.single_scattering_albedo_;
  ArrayPolicy& assymetry_parameter = accumulated_radiator_states.assymetry_parameter_;
  for (std::size_t i = 0; i < number_of_columns; i++)
  {
    // delta eddington parameters
    auto& gamma1 = parameters[i][0];
    auto& gamma2 = parameters[i][1];
    auto& gamma3 = parameters[i][2];
    auto& gamma4 = parameters[i][3];
    auto& mu = parameters[i][4];

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
  }

  template<typename T, typename GridPolicy, typename ProfilePolicy, typename RadiatorStatePolicy>
  void AssembleTridiagonalSystem(
      const std::vector<T>& solar_zenith_angles,
      const std::map<std::string, GridPolicy>& grids,
      const std::map<std::string, ProfilePolicy>& profiles,
      const std::vector<T> solution_parameters,
      const TridiagonalMatrix<T>& coeffcient_matrix,
      const std::vector<T>& coeffcient_vector)
  {
    // get system size
    std::size_t system_size = coeffcient_matrix.main_diagonal_.size() / 2;

    // coeffcient matrix diagonals
    std::vector<T>& upper_diagonal = coeffcient_matrix.upper_diagonal_;
    std::vector<T>& main_diagonal = coeffcient_matrix.main_diagonal_;
    std::vector<T>& lower_diagonal = coeffcient_matrix.lower_diagonal;

    {
      // first row
      lower_diagonal[0] = 0;
      main_diagonal[0] = solution_parameters[0][0];
      upper_diagonal[0] = -solution_parameters[1][0];

      // even rows
      for (std::size_t n = 2; n < 2 * system_size / 2; n += 2)
      {
        const auto& e1 = solution_parameters[n][0];
        const auto& e2 = solution_parameters[n][1];
        const auto& e3 = solution_parameters[n][2];
        const auto& e4 = solution_parameters[n][3];
        upper_diagonal[n] = e2[n + 1] * e1[n] - e3[n] * e4[n + 1];
        main_diagonal[n] = e2[n] * e2[n + 1] - e4[n] * e4[n + 1];
        lower_diagonal[n] = e1[n] * e4[n + 1] - e2[n + 1] * e3[n + 1];
      }

      // odd rows
      for (std::size_t n = 1; n < 2 * system_size / 2 - 1; n += 2)
      {
        const auto& e1 = solution_parameters[n][0];
        const auto& e2 = solution_parameters[n][1];
        const auto& e3 = solution_parameters[n][2];
        const auto& e4 = solution_parameters[n][3];
        upper_diagonal[n] = e2[n] * e3[n] - e4[n] * e1[n];
        main_diagonal[n] = e1[n] * e1[n + 1] - e3[n] * e3[n + 1];
        lower_diagonal[n] = e3[n] * e4[n + 1] - e1[n + 1] * e2[n + 1];
      }

      // last row
      const auto& e1 = solution_parameters[n][0];
      const auto& e2 = solution_parameters[n][1];
      const auto& e3 = solution_parameters[n][2];
      const auto& e4 = solution_parameters[n][3];
      std::size_t N = 2 * system_size - 1;
      upper_diagonal[N] = e1[N];
    }
  }

  /// @brief Initilize the variables required for the delta eddington approximation
  /// @param solar zenith angles Solar zenith angles for each column [radians]
  /// @param grids Grids available for the radiative transfer calculation
  /// @param profiles Profiles available for the radiative transfer calculation
  /// @param radiation_field The calculated radiation field
  ///
  /// Delta scaling of the inputs, compuation of reflectivity, diffuse flux and
  /// delta eddington coeffcients.
  template<
      typename T,
      typename GridPolicy,
      typename ProfilePolicy,
      typename RadiatorStatePolicytypename,
      typename RadiationFieldPolicy>
  void ComputeRadiationField(
      const std::vector<T>& solar_zenith_angles,
      const std::map<std::string, GridPolicy>& grids,
      const std::map<std::string, ProfilePolicy>& profiles,
      const Array2D<T> solution_parameters,
      const TridiagonalMatrix<T>& coeffcient_matrix,
      const std::vector<T>& coeffcient_vector,
      const RadiationFieldPolicy& radiation_field)
  {
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
    Array2D<T> simulation_parameters;

    // tridiagonal system variables
    TridiagonalMatrix<T> coeffcient_matrix;
    std::vector<T> coeffcient_vector;

    tuvx::InitializeVariables<T, GridPolicy, ProfilePolicy, RadiatorStatePolicy, RadiationFieldPolicy>(
        solar_zenith_angles, grids, profiles, accumulated_radiator_state);

    ApproximationFunction(accumulated_radiator_state, solar_zenith_angles, simulation_parameters);

    tuvx::AssembleTridiagonalSystem<T, GridPolicy, ProfilePolicy, RadiatorStatePolicy, RadiationFieldPolicy>(
        solar_zenith_angles, grids, profiles, solution_parameters, coeffcient_matrix, coeffcient_vector);

    tuvx::Solve<T>(coeffcient_matrix, coeffcient_vector);

    tuvx::ComputeRadiationField<T, GridPolicy, ProfilePolicy, RadiatorStatePolicy>(
        solar_zenith_angles, grids, profiles, solution_parameters, coeffcient_matrix, coeffcient_vector, radiation_field);
  }

}  // namespace tuvx
