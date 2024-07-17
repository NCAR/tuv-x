
// Copyright (C) 2023-2024 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#include "tuvx/linear_algebra/linear_algebra.hpp"
#include "tuvx/radiative_transfer/radiator.hpp"

#include <tuvx/radiative_transfer/solvers/delta_eddington.hpp>

#include <cmath>
#include <functional>

namespace tuvx
{

  template<
      typename T,
      typename GridPolicy,
      typename ProfilePolicy,
      typename RadiatorStatePolicy,
      typename RadiationFieldPolicy,
      typename SourceFunctionPolicy,
      typename ArrayPolicy>
  inline void InitializeVariables(
      const std::vector<T>& solar_zenith_angles,
      const std::map<std::string, GridPolicy>& grids,
      const std::map<std::string, ProfilePolicy>& profiles,
      const std::map<std::string, ArrayPolicy>& solver_parameters,
      const std::map<std::string, ArrayPolicy>& solution_parameters,
      std::map<std::string, ArrayPolicy>& source_functions,
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
    ArrayPolicy& tau = accumulated_radiator_states.optical_depth_;
    ArrayPolicy& omega = accumulated_radiator_states.single_scattering_albedo_;
    ArrayPolicy& g = accumulated_radiator_states.assymetry_parameter_;
  }

  template<typename T>
  void ScaleVariables(std::map<std::string, std::vector<T>>& solver_parameters)
  {
    // solver parameters (environment related variables)
    std::vector<T>& solar_zenith_angles = solver_parameters.at("Solar Zenith Angles");
    std::vector<T>& tau = solver_parameters.at("Optical Depth");
    std::vector<T>& g = solver_parameters.at("Assymetry Parameter");
    std::vector<T>& omega = solver_parameters.at("Single Scattering Albedo");
    std::size_t number_of_columns = solar_zenith_angles.size();

    // Delta scaling
    T f;
    for (std::size_t i = 0; i < number_of_columns; i++)
    {
      f = omega[i] * omega[i];
      g[i] = (g[i] - f) / (1 - f);
      omega[i] = (1 - f) * omega[i] / (1 - omega[i] * f);
      tau[i] = (1 - omega[i] * f) * tau[i];
    }

    // TODO slant optical depth computation
    for (auto& tau_n : tau)
    {
      tau_n = tau_n;
    }
  }

  template<typename T>
  void BuildSourceFunctions(
      std::map<std::string, std::vector<T>> solver_parameters,
      std::map<std::string, std::vector<T>> solution_parameters,
      std::map<std::string, std::function<T(T)>> source_functions)
  {
    // Source terms (C1 and C2 from the paper)
    auto& C_upwelling = source_functions.at("C_upwelling");
    auto& C_downwelling = source_functions.at("C_downwelling");

    // solver parameters (environment related variables)
    std::vector<T>& solar_zenith_angles = solver_parameters.at("Solar Zenith Angles");
    std::vector<T>& tau = solver_parameters.at("Optical Depth");
    std::vector<T>& omega = solver_parameters.at("Assymetry Parameter");
    std::size_t number_of_columns = solar_zenith_angles.size();

    // parameters used to compute the soulution
    std::vector<T>& lambda = solution_parameters.at("lambda");
    std::vector<T>& gamma1 = solution_parameters.at("gamma1");
    std::vector<T>& gamma2 = solution_parameters.at("gamma2");
    std::vector<T>& gamma3 = solution_parameters.at("gamma3");
    std::vector<T>& gamma4 = solution_parameters.at("gamma4");
    std::vector<T>& mu = solution_parameters.at("mu");

    // solution parameters
    auto& S_sfc_i = solution_parameters.at("Infrared Source Flux");
    auto& S_sfc_s = solution_parameters.at("Solar Source Flux");
    auto& R_sfc = solution_parameters.at("source flux");

    // temporary variables
    T tau_cumulative = 0;
    T exponential_term, denominator_term, mu_0;

    // source terms (C equations from 16, 290; eqns 23, 24)
    for (std::size_t i = 0; i < number_of_columns; i++)
    {
      mu_0 = std::acos(solar_zenith_angles[i]);
      denominator_term = (lambda * lambda - 1 / (mu_0 * mu_0));
      tau_cumulative += tau[i];

      S_sfc_i[i] = R_sfc * mu_0 * std::exp(-tau_cumulative / mu_0);
      S_sfc_s[i] = M_PI * R_sfc;

      C_downwelling[i] = [&i, &mu_0](T tau) -> T
      {
        T exponential_term = omega * M_PI * R_sfc * std::exp(-(tau_cumulative - tau) / mu_0);
        return exponential_term * (((gamma1[i] + 1) / mu_0) * gamma4[i] + gamma2[i] * gamma3[i]);
      };
      C_upwelling[i] = [&i, &mu_0](T tau) -> T
      {
        T exponential_term = omega * M_PI * R_sfc * std::exp(-(tau_cumulative - tau) / mu_0);
        return exponential_term * (((gamma1[i] - 1) / mu_0) * gamma3[i] + gamma4[i] * gamma2[i]);
      };
    }
  }
}  // namespace tuvx
