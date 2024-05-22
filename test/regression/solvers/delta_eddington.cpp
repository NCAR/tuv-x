// Copyright (C) 2023-2024 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
//
// Tests for tuvx::DeltaEddington.
#include "delta_eddington.hpp"
#include <iostream>

#define ASSERT(x) if (!(x)) { \
    std::cerr << "Assertion failed at " << __FILE__ << ":" << __LINE__ << " : "<< #x << std::endl; \
    exit(EXIT_FAILURE); \
}

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

// Checks certian input values to ensure the data is properly transferred
// from the Fortran side to the C++ side.
// If the testing conditions change in the future, this function will need to be updated.
void CheckInputs(const std::vector<double>& solar_zenith_angles,
                 const std::vector<double>& earth_sun_distances,
                 const std::map<std::string, tuvx::Grid<tuvx::Array2D<double>>>& grids,
                 const std::map<std::string, tuvx::Profile<tuvx::Array2D<double>>>& profiles,
                 const tuvx::RadiatorState<tuvx::Array3D<double>>& accumulated_radiator_states)
{
  ASSERT(solar_zenith_angles.size() == 2);
  ASSERT(earth_sun_distances.size() == 2);
  ASSERT(grids.size() == 2);
  ASSERT(grids.at("altitude [m]").NumberOfColumns() == 2);
  ASSERT(grids.at("altitude [m]").NumberOfSections() == 120);
  ASSERT(!grids.at("altitude [m]").IsConstant());
  // heights have been converted to meters from km
  ASSERT(grids.at("altitude [m]").edges_(0,0) == 0.0); 
  ASSERT(grids.at("altitude [m]").edges_(42,0) == 42000.0);
  ASSERT(grids.at("altitude [m]").edges_(120,0) == 120000.0); 
  ASSERT(grids.at("altitude [m]").edges_(0,1) == 0.0); 
  ASSERT(grids.at("altitude [m]").edges_(42,1) == 42000.0);
  ASSERT(grids.at("altitude [m]").edges_(120,1) == 120000.0); 
  ASSERT(grids.at("altitude [m]").mid_points_(0,0) == 500.0);
  ASSERT(grids.at("altitude [m]").mid_points_(41,0) == 41500.0);
  ASSERT(grids.at("altitude [m]").mid_points_(119,0) == 119500.0);
  ASSERT(grids.at("altitude [m]").mid_points_(0,1) == 500.0);
  ASSERT(grids.at("altitude [m]").mid_points_(41,1) == 41500.0);
  ASSERT(grids.at("altitude [m]").mid_points_(119,1) == 119500.0);
  ASSERT(grids.at("wavelength [m]").NumberOfColumns() == 1);
  ASSERT(grids.at("wavelength [m]").NumberOfSections() == 156);
  ASSERT(grids.at("wavelength [m]").IsConstant());
  // wavelengths have been converted to meters from nm
  ASSERT(grids.at("wavelength [m]").edges_(0,0) == 120.0*1.0e-9);
  ASSERT(grids.at("wavelength [m]").edges_(77,0) == 311.5*1.0e-9);
  ASSERT(grids.at("wavelength [m]").edges_(156,0) == 735.0*1.0e-9);
  ASSERT(grids.at("wavelength [m]").mid_points_(0,0) > 120.69*1.0e-9);
  ASSERT(grids.at("wavelength [m]").mid_points_(0,0) < 120.71*1.0e-9);
  ASSERT(grids.at("wavelength [m]").mid_points_(77,0) > 311.9*1.0e-9);
  ASSERT(grids.at("wavelength [m]").mid_points_(77,0) < 312.1*1.0e-9);
  ASSERT(grids.at("wavelength [m]").mid_points_(155,0) > 729.9*1.0e-9);
  ASSERT(grids.at("wavelength [m]").mid_points_(155,0) < 730.1*1.0e-9);
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
  CheckInputs(solar_zenith_angles, earth_sun_distances, grids, profiles, accumulated_radiator_states);
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
