
// Copyright (C) 2023-2024 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

namespace tuvx
{

  template<
      typename T,
      typename GridPolicy,
      typename RadiatorStatePolicy,
      typename RadiationFieldPolicy,
      typename SourceFunctionPolicy,
      typename ArrayPolicy>
  inline void InitializeVariables(
      const std::vector<T>& solar_zenith_angles,
      const std::map<std::string, GridPolicy>& grids,
      const RadiatorStatePolicy& accumulated_radiator_states,
      std::map<std::string, ArrayPolicy>& solver_variables,
      std::map<std::string, SourceFunctionPolicy>& source_functions)
  {
    // determine number of layers
    const std::size_t number_of_columns = solar_zenith_angles.size();
    const auto& vertical_grid = grids.at("altitude [m]");
    const auto& wavelength_grid = grids.at("wavelength [m]");

    // Check for consistency between the grids and profiles.
    assert(vertical_grid.NumberOfColumns() == number_of_columns);
    assert(wavelength_grid.NumberOfColumns() == 1);

    // Scale variables for the Delta-Eddington approximation
    ScaleVariables(grids, solar_zenith_angles, solver_variables);

    // Generate functions that define the source variables (C functions from the paper)
    BuildSourceFunctions(grids, solar_zenith_angles, solver_variables, source_functions);
  }

  template<typename T, typename ArrayPolicy, typename GridPolicy>
  inline void ScaleVariables(
      GridPolicy grids,
      std::vector<T> solar_zenith_angles,
      std::map<std::string, std::vector<T>>& solver_variables)
  {
    // solver parameters (environment related variables)
    std::vector<T>& tau = solver_variables.at("Optical Depth");
    std::vector<T>& g = solver_variables.at("Assymetry Parameter");
    std::vector<T>& omega = solver_variables.at("Single Scattering Albedo");

    // grid and dimensions
    const auto& vertical_grid = grids.at("altitude [m]");
    const auto& wavelength_grid = grids.at("wavelength [m]");

    // dimensions
    const std::size_t number_of_columns = solar_zenith_angles.size();
    const std::size_t number_of_wavelengths = wavelength_grid.NumberOfColumns();
    const std::size_t number_of_layers = vertical_grid.NumberOfColumns();

    // Delta scaling
    T f;
    for (std::size_t i = 0; i < number_of_layers; i++)
    {
      for (std::size_t j = 0; j < number_of_wavelengths; j++)
      {
        for (std::size_t k = 0; k < number_of_columns; k++)
        {
          f = omega[i][j][j] * omega[i][j][k];
          g[i][j][k] = (g[i][j][k] - f) / (1 - f);
          omega[i][j][k] = (1 - f) * omega[i][j][k] / (1 - omega[i][j][k] * f);
          tau[i][j][k] = (1 - omega[i][j][k] * f) * tau[i][j][k];
        }
      }
    }

    // TODO slant optical depth computation
    // for (auto& tau_n : tau)
    //{
    //  tau_n = tau_n;
    //}
  }

  template<typename T, typename GridPolicy, typename ArrayPolicy, typename SourceFunctionPolicy>
  inline void BuildSourceFunctions(
      GridPolicy grids,
      std::vector<T> solar_zenith_angles,
      std::map<std::string, std::vector<T>> solver_variables,
      std::map<std::string, SourceFunctionPolicy> source_functions)
  {

    // grid and dimensions
    const auto& vertical_grid = grids.at("altitude [m]");
    const auto& wavelength_grid = grids.at("wavelength [m]");

    // dimensions
    const std::size_t& number_of_columns = solar_zenith_angles.size();
    const std::size_t& number_of_wavelengths = wavelength_grid.NumberOfColumns();
    const std::size_t& number_of_layers = vertical_grid.NumberOfColumns();

    // solver parameters (environment related variables)
    const auto& tau = solver_variables.at("Optical Depth");
    const auto& omega = solver_variables.at("Assymetry Parameter");

    // parameters used to compute the soulution
    const auto& lambda = solver_variables.at("lambda");
    const auto& gamma1 = solver_variables.at("gamma1");
    const auto& gamma2 = solver_variables.at("gamma2");
    const auto& gamma3 = solver_variables.at("gamma3");
    const auto& gamma4 = solver_variables.at("gamma4");
    const auto& mu = solver_variables.at("mu");

    // solution parameters
    auto& S_sfc_i = solver_variables.at("Infrared Source Flux");
    auto& S_sfc_s = solver_variables.at("Solar Source Flux");
    auto& R_sfc = solver_variables.at("source flux");

    // temporary variables
    T tau_cumulative = 0;
    T exponential_term, denominator_term, mu_0;

    // Source terms (C1 and C2 from the paper)
    auto& C_upwelling = source_functions.at("C_upwelling");
    auto& C_downwelling = source_functions.at("C_downwelling");

    // source terms (C equations from 16, 290; eqns 23, 24)
    for (std::size_t i = 0; i < number_of_layers; i++)
    {
      for (std::size_t j = 0; j < number_of_wavelengths; j++)
      {
        for (std::size_t k = 0; k < number_of_columns; k++)
        {

          mu_0 = std::acos(solar_zenith_angles[i][j][k]);
          denominator_term = (lambda[i][j][k] * lambda[i][j][k] - 1 / (mu_0 * mu_0));
          tau_cumulative += tau[i][j][k];

          S_sfc_i[i][j][k] = R_sfc[i][j][k] * mu_0 * std::exp(-tau_cumulative / mu_0);
          S_sfc_s[i][j][k] = M_PI * R_sfc[i][j][k];

          // Source function defined for each grid point
          C_downwelling[i][j][k] = [&i, &j, &k, &mu_0](T tau) -> T
          {
            T exponential_term = omega[i][j][k] * M_PI * R_sfc[i][j][k] * std::exp(-(tau_cumulative - tau) / mu_0);
            return exponential_term * (((gamma1[i][j][k] + 1) / mu_0) * gamma4[i][j][k] + gamma2[i][j][k] * gamma3[i][j][k]);
          };

          C_upwelling[i][j][k] = [&i, &j, &k, &mu_0](T tau) -> T
          {
            T exponential_term = omega[i][j][k] * M_PI * R_sfc[i][j][k] * std::exp(-(tau_cumulative - tau) / mu_0);
            return exponential_term * (((gamma1[i][j][k] - 1) / mu_0) * gamma3[i][j][k] + gamma4[i][j][j] * gamma2[i][j][k]);
          };

        }
      }
    }

  }

}  // namespace tuvx
