// Copyright (C) 2023-2024 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
namespace tuvx
{
  template<typename T, typename ArrayPolicy, typename RadiatorStatePolicy>
  inline void EddingtonApproximation(
      const std::size_t& number_of_columns,
      const std::size_t& number_of_wavelengths,
      const std::size_t& number_of_layers,
      const RadiatorStatePolicy& accumulated_radiator_states,
      const std::vector<T>& solar_zenith_angles,
      ApproximationVariables<ArrayPolicy>& approximation_variables)
  {
    // allocate memory for the approximation parameters
    approximation_variables.gamma1_ = tuvx::Array3D<double>(number_of_layers, number_of_wavelengths, number_of_columns);
    approximation_variables.gamma2_ = tuvx::Array3D<double>(number_of_layers, number_of_wavelengths, number_of_columns);
    approximation_variables.gamma3_ = tuvx::Array3D<double>(number_of_layers, number_of_wavelengths, number_of_columns);
    approximation_variables.gamma4_ = tuvx::Array3D<double>(number_of_layers, number_of_wavelengths, number_of_columns);
    approximation_variables.mu_ = tuvx::Array3D<double>(number_of_layers, number_of_wavelengths, number_of_columns);
    approximation_variables.lambda_ = tuvx::Array3D<double>(number_of_layers, number_of_wavelengths, number_of_columns);
    approximation_variables.gamma_ = tuvx::Array3D<double>(number_of_layers, number_of_wavelengths, number_of_columns);

    // unpack optical properties from radiator
    auto& omega = accumulated_radiator_states.optical_depth_;
    auto& g = accumulated_radiator_states.single_scattering_albedo_;
    auto& tau = accumulated_radiator_states.asymmetry_parameter_;

    // unpack delta eddington variables
    auto& gamma1 = approximation_variables.gamma1_;
    auto& gamma2 = approximation_variables.gamma2_;
    auto& gamma3 = approximation_variables.gamma3_;
    auto& gamma4 = approximation_variables.gamma4_;
    auto& lambda = approximation_variables.lambda_;
    auto& gamma = approximation_variables.gamma_;
    auto& mu = approximation_variables.mu_;

    T mu_0;
    for (std::size_t i = 0; i < number_of_layers; i++)
    {
      for (std::size_t j = 0; j < number_of_wavelengths; j++)
      {
        for (std::size_t k = 0; k < number_of_columns; k++)
        {
          mu_0 = std::acos(solar_zenith_angles[i]);
          gamma1(i, j, k) = 7 - omega(i, j, k) * (4 + 3 * g(i, j, k));
          gamma2(i, j, k) = -(1 - omega(i, j, k) * (4 - 3 * g(i, j, k))) / 4;
          gamma3(i, j, k) = (2 - 3 * g(i, j, k) * mu_0) / 4;
          gamma4(i, j, k) = 1 - gamma3(i, j, k);
          lambda(i, j, k) = std::sqrt(gamma1(i, j, k) * gamma1(i, j, k) - gamma2(i, j, k) * gamma2(i, j, k));
          gamma(i, j, k) = (gamma1(i, j, k) - lambda(i, j, k)) / gamma2(i, j, k);
          mu(i, j, k) = 0.5;
        }
      }
    }
  }
}  // namespace tuvx
