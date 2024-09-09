
// Copyright (C) 2023-2024 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0

namespace tuvx
{

  template<typename T, typename GridPolicy, typename RadiatorStatePolicy, typename SourceFunctionPolicy>
  inline void InitializeSolver(
      const std::vector<T>& solar_zenith_angles,
      const std::map<std::string, GridPolicy>& grids,
      RadiatorStatePolicy& accumulated_radiator_states,
      std::map<std::string, SourceFunctionPolicy>& source_functions)
  {
    // determine number of layers
    const std::size_t number_of_columns = solar_zenith_angles.size();
    const auto& vertical_grid = grids.at("altitude [m]");
    const auto& wavelength_grid = grids.at("wavelength [m]");

    // Check for consistency between the grids and profiles.
    assert(vertical_grid.NumberOfColumns() == number_of_columns);

    // Scale variables for the Delta-Eddington approximation
    ScaleVariables(grids, solar_zenith_angles, accumulated_radiator_states);

    // Generate functions that define the source variables (C functions from the paper)
    // BuildSourceFunctions(grids, solar_zenith_angles, solver_variables, source_functions);
  }

  template<typename T, typename GridPolicy, typename RadiatorStatePolicy>
  inline void ScaleVariables(
      const std::map<std::string, GridPolicy>& grids,
      const std::vector<T>& solar_zenith_angles,
      RadiatorStatePolicy& radiator_state)
  {
    // solver parameters (environment related variables)
    auto& tau = radiator_state.optical_depth_;
    auto& g = radiator_state.asymmetry_parameter_;
    auto& omega = radiator_state.single_scattering_albedo_;

    // grids
    const auto& vertical_grid = grids.at("altitude [m]");
    const auto& wavelength_grid = grids.at("wavelength [m]");

    // dimensions
    const std::size_t& number_of_columns = solar_zenith_angles.size();
    const std::size_t& number_of_wavelengths = wavelength_grid.NumberOfColumns();
    const std::size_t& number_of_layers = vertical_grid.NumberOfColumns();

    // Delta scaling
    T f;
    for (std::size_t i = 0; i < number_of_layers; i++)
    {
      for (std::size_t j = 0; j < number_of_wavelengths; j++)
      {
        for (std::size_t k = 0; k < number_of_columns; k++)
        {
          f = omega(i, j, k) * omega(i, j, k);
          g(i, j, k) = (g(i, j, k) - f) / (1 - f);
          omega(i, j, k) = (1 - f) * omega(i, j, k) / (1 - omega(i, j, k) * f);
          tau(i, j, k) = (1 - omega(i, j, k) * f) * tau(i, j, k);
        }
      }
    }

    // TODO slant optical depth computation
    // for (auto& tau_n : tau)
    //{
    //  tau_n = tau_n;
    //}
  }

  template<typename T, typename ArrayPolicy, typename RadiatorStatePolicy, typename SourceFunctionPolicy>
  inline void BuildSourceFunctions(
      const std::size_t& number_of_columns,
      const std::size_t& number_of_wavelengths,
      const std::size_t& number_of_layers,
      const std::vector<T>& solar_zenith_angles,
      const RadiatorStatePolicy& accumulated_radiator_states,
      const ApproximationVariables<ArrayPolicy> approximation_variables,
      const ArrayPolicy& surface_reflectivity,
      std::map<std::string, SourceFunctionPolicy> source_functions)
  {
    // solver parameters (environment related variables)
    const auto& tau = accumulated_radiator_states.optical_depth_;
    const auto& omega = accumulated_radiator_states.asymmetry_parameter_;

    // parameters used to compute the soulution
    const auto& mu = approximation_variables.mu_;
    const auto& lambda = approximation_variables.lambda_;
    const auto& gamma1 = approximation_variables.gamma1_;
    const auto& gamma2 = approximation_variables.gamma2_;
    const auto& gamma3 = approximation_variables.gamma3_;
    const auto& gamma4 = approximation_variables.gamma4_;
    const auto& R_sfc = surface_reflectivity;

    // temporary variables
    T tau_cumulative = 0;
    T exponential_term, denominator_term, mu_0;

    // Source terms (C1 and C2 from the paper)
    auto& C_upwelling = source_functions.at("C_upwelling");
    auto& C_downwelling = source_functions.at("C_downwelling");

    for (std::size_t i = 0; i < number_of_columns; i++)
    {
      mu_0 = std::acos(solar_zenith_angles[i]);
      for (std::size_t j = 0; j < number_of_wavelengths; j++)
      {
        for (std::size_t k = 0; k < number_of_layers; k++)
        {
          denominator_term = (lambda(i, j, k) * lambda(i, j, k) - 1 / (mu_0 * mu_0));
          tau_cumulative += tau(i, j, k);

          // Source function defined for each grid point
          C_downwelling(i, j, k) =
              [&i, &j, &k, &mu_0, &R_sfc, &gamma1, &gamma2, &gamma3, &gamma4, &lambda, &mu, &omega, &tau_cumulative](
                  T tau) -> T
          {
            T exponential_term = omega(i, j, k) * M_PI * R_sfc(i, j, k) * std::exp(-(tau_cumulative - tau) / mu_0);
            return exponential_term * (((gamma1(i, j, k) + 1) / mu_0) * gamma4(i, j, k) + gamma2(i, j, k) * gamma3(i, j, k));
          };

          C_upwelling(i, j, k) =
              [&i, &j, &k, &mu_0, &R_sfc, &gamma1, &gamma2, &gamma3, &gamma4, &lambda, &mu, &omega, &tau_cumulative](
                  T tau) -> T
          {
            T exponential_term = omega(i, j, k) * M_PI * R_sfc(i, j, k) * std::exp(-(tau_cumulative - tau) / mu_0);
            return exponential_term * (((gamma1(i, j, k) - 1) / mu_0) * gamma3(i, j, k) + gamma4(i, j, k) * gamma2(i, j, k));
          };
        }
      }
    }
  }

}  // namespace tuvx
