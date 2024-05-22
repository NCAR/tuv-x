// Copyright (C) 2023-2024 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
//
// Tests for tuvx::DeltaEddington.
#include "delta_eddington.hpp"

double* CopyVector(const std::vector<double>& vec)
{
  double* arr = new double[vec.size()];
  for (int i = 0; i < vec.size(); i++) {
    arr[i] = vec[i];
  }
  return arr;
}

template<typename GridPolicy>
GridPolicy CreateGrid(std::string units, const size_t n_columns, const size_t n_sections, const double* mid_points, const double* edges)
{
  GridPolicy grid(units, n_columns, n_sections);
  size_t index = 0;
  for (auto &elem : grid.mid_points_) {
    elem = mid_points[index++];
  }
  index = 0;
  for (auto &elem : grid.edges_) {
    elem = edges[index++];
  }
  return grid;
}

template<typename GridPolicy>
GridPolicy CreateFixedGrid(std::string units, const size_t sections, const double* mid_points, const double* edges)
{
  GridPolicy grid(units, sections);
  size_t index = 0;
  for (auto &elem : grid.mid_points_) {
    elem = mid_points[index++];
  }
  index = 0;
  for (auto &elem : grid.edges_) {
    elem = edges[index++];
  }
  return grid;
}

SolverOutput RunDeltaEddingtonSolver(const SolverInput input)
{
  using GridPolicy = tuvx::Grid<tuvx::Array2D<double>>;
  std::vector<double> solar_zenith_angles(input.solar_zenith_angles_, input.solar_zenith_angles_ + input.n_columns_);
  std::vector<double> earth_sun_distances(input.earth_sun_distances_, input.earth_sun_distances_ + input.n_columns_);
  std::map<std::string, GridPolicy> grids;
  grids["altitude [m]"]   = CreateGrid<GridPolicy>("m", input.n_columns_, input.n_levels_, input.altitude_mid_points_, input.altitude_edges_);
  grids["wavelength [m]"] = CreateFixedGrid<GridPolicy>("m", input.n_wavelengths_, input.wavelength_mid_points_, input.wavelength_edges_);
  std::map<std::string, tuvx::Profile<tuvx::Array2D<double>>> profiles;
  tuvx::RadiatorState<tuvx::Array3D<double>> accumulated_radiator_states(input.n_columns_,
                                                                         grids["altitude [m]"],
                                                                         grids["wavelength [m]"]);
  tuvx::RadiationField<tuvx::RadiationFieldComponents<tuvx::Array3D<double>>> radiation_field(input.n_columns_,
                                                                                              grids["altitude [m]"],
                                                                                              grids["wavelength [m]"]);
  tuvx::DeltaEddington solver;
  solver.Solve(solar_zenith_angles, grids, profiles, accumulated_radiator_states, radiation_field);
  SolverOutput output;
  output.n_wavelengths_ = input.n_wavelengths_;
  output.n_levels_ = input.n_levels_;
  output.n_columns_ = input.n_columns_;
  output.flux_direct_ = CopyVector(radiation_field.actinic_flux_.direct_.AsVector());
  output.flux_up_ = CopyVector(radiation_field.actinic_flux_.upwelling_.AsVector());
  output.flux_down_ = CopyVector(radiation_field.actinic_flux_.downwelling_.AsVector());
  output.irrad_direct_ = CopyVector(radiation_field.spectral_irradiance_.direct_.AsVector());
  output.irrad_up_ = CopyVector(radiation_field.spectral_irradiance_.upwelling_.AsVector());
  output.irrad_down_ = CopyVector(radiation_field.spectral_irradiance_.downwelling_.AsVector());
  return output;
}

void FreeOutput(SolverOutput output)
{
  delete[] output.flux_direct_;
  delete[] output.flux_up_;
  delete[] output.flux_down_;
  delete[] output.irrad_direct_;
  delete[] output.irrad_up_;
  delete[] output.irrad_down_;
}
