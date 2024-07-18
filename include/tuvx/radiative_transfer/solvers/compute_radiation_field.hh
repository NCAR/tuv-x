

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
      typename GridPolicy,
      typename ProfilePolicy,
      typename RadiatorStatePolicytypename,
      typename RadiationFieldPolicy>
  void ComputeRadiationField(
      const std::vector<T>& solar_zenith_angles,
      const std::map<std::string, GridPolicy>& grids,
      const std::map<std::string, ProfilePolicy>& profiles,
      const std::map<std::string, std::vector<T>> solution_parameters,
      RadiationFieldPolicy& radiation_field)
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

}  // namespace tuvx
