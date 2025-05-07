// Copyright (C) 2023-2025 University Corporation for Atmospheric Research-National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <tuvx/radiative_transfer/solvers/delta_eddington.hpp>

#ifdef __cplusplus
extern "C"
{
#endif

  struct SolverInput
  {
    int n_wavelengths_;
    int n_levels_;
    int n_columns_;
    double* solar_zenith_angles_;        // [columns]
    double* earth_sun_distances_;        // [columns]
    double* altitude_mid_points_;        // [levels][columns]
    double* altitude_edges_;             // [levels + 1][columns]
    double* wavelength_mid_points_;      // [wavelengths]
    double* wavelength_edges_;           // [wavelengths + 1]
    double* optical_depths_;             // [wavelengths][levels][columns]
    double* single_scattering_albedos_;  // [wavelengths][levels][columns]
    double* asymmetry_parameters_;       // [wavelengths][levels][columns]
  };

  struct SolverOutput
  {
    int n_wavelengths_;
    int n_levels_;
    int n_columns_;
    double* flux_direct_;   // [wavelengths][levels+1][columns]
    double* flux_up_;       // [wavelengths][levels+1][columns]
    double* flux_down_;     // [wavelengths][levels+1][columns]
    double* irrad_direct_;  // [wavelengths][levels+1][columns]
    double* irrad_up_;      // [wavelengths][levels+1][columns]
    double* irrad_down_;    // [wavelengths][levels+1][columns]
  };

  SolverOutput RunDeltaEddingtonSolver(const SolverInput input);

  void FreeOutput(SolverOutput output);

#ifdef __cplusplus
}  // extern "C"
#endif