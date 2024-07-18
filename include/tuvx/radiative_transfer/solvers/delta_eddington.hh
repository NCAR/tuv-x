// Copyright (C) 2023-2024 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#include "tuvx/linear_algebra/linear_algebra.hpp"
#include "tuvx/radiative_transfer/radiator.hpp"

#include <tuvx/radiative_transfer/solvers/delta_eddington.hpp>

#include <cmath>

namespace tuvx
{
  template<typename T, typename ArrayPolicy, typename RadiatorStatePolicy>
  void DeltaEddingtonApproximation(
      const RadiatorStatePolicy& accumulated_radiator_states,
      const std::map<std::string, std::vector<T>> solution_parameters,
      const std::vector<T> solar_zenith_angles)
  {
    const std::size_t number_of_columns = solar_zenith_angles.size();
    ArrayPolicy& omega = accumulated_radiator_states.optical_depth_;
    ArrayPolicy& g = accumulated_radiator_states.single_scattering_albedo_;
    ArrayPolicy& tau = accumulated_radiator_states.assymetry_parameter_;

    // delta eddington parameters
    std::vector<T>& gamma1 = solution_parameters.at("gamma1");
    std::vector<T>& gamma2 = solution_parameters.at("gamma2");
    std::vector<T>& gamma3 = solution_parameters.at("gamma3");
    std::vector<T>& gamma4 = solution_parameters.at("gamma4");
    std::vector<T>& mu = solution_parameters.at("mu");
    std::vector<T>& lambda = solution_parameters.at("lambda");
    std::vector<T>& gamma = solution_parameters.at("Gamma");

    // simulation parameters
    T mu_0;
    for (std::size_t i = 0; i < number_of_columns; i++)
    {
      // compute delta eddington parameters
      mu_0 = std::acos(solar_zenith_angles[i]);
      gamma1[i] = 7 - omega[i] * (4 + 3 * g[i]);
      gamma2[i] = -(1 - omega[i] * (4 - 3 * g[i])) / 4;
      gamma3[i] = (2 - 3 * g[i] * mu_0) / 4;
      gamma4[i] = 1 - gamma3[i];
      lambda[i] = std::sqrt(gamma1[i] * gamma1[i] - gamma2[i] * gamma2[i]);
      gamma[i] = (gamma1[i] - lambda[i]) / gamma2[i];
      mu = (T)0.5;
    }
  }
}  // namespace tuvx
