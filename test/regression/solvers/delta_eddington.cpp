// Copyright (C) 2023-2024 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
//
// Tests for tuvx::DeltaEddington.
#include <gtest/gtest.h>

#include "delta_eddington.hpp"
#include <iostream>

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
  // Solar zenith angles have been converted to radians from degrees
  ASSERT_EQ(solar_zenith_angles.size(), 2);
  ASSERT_NEAR(solar_zenith_angles.at(0), 1.8294 * 3.14159265 / 180.0, 1.0e-6);
  ASSERT_NEAR(solar_zenith_angles.at(1), 28.199 * 3.14159265 / 180.0, 1.0e-6);

  // Earth-Sun distances are in AU
  ASSERT_EQ(earth_sun_distances.size(), 2);
  ASSERT_NEAR(earth_sun_distances.at(0), 0.9961, 1.0e-6);
  ASSERT_NEAR(earth_sun_distances.at(1), 0.9962, 1.0e-6);

  // heights have been converted to meters from km
  ASSERT_EQ(grids.size(), 2);
  ASSERT_EQ(grids.at("altitude [m]").NumberOfColumns(), 2);
  ASSERT_EQ(grids.at("altitude [m]").NumberOfSections(), 120);
  ASSERT_FALSE(grids.at("altitude [m]").IsConstant());
  ASSERT_EQ(grids.at("altitude [m]").edges_(0,0), 0.0); 
  ASSERT_EQ(grids.at("altitude [m]").edges_(42,0), 42000.0);
  ASSERT_EQ(grids.at("altitude [m]").edges_(120,0), 120000.0); 
  ASSERT_EQ(grids.at("altitude [m]").edges_(0,1), 0.0); 
  ASSERT_EQ(grids.at("altitude [m]").edges_(42,1), 42000.0);
  ASSERT_EQ(grids.at("altitude [m]").edges_(120,1), 120000.0); 
  ASSERT_EQ(grids.at("altitude [m]").mid_points_(0,0), 500.0);
  ASSERT_EQ(grids.at("altitude [m]").mid_points_(41,0), 41500.0);
  ASSERT_EQ(grids.at("altitude [m]").mid_points_(119,0), 119500.0);
  ASSERT_EQ(grids.at("altitude [m]").mid_points_(0,1), 500.0);
  ASSERT_EQ(grids.at("altitude [m]").mid_points_(41,1), 41500.0);
  ASSERT_EQ(grids.at("altitude [m]").mid_points_(119,1), 119500.0);

  // wavelengths have been converted to meters from nm
  ASSERT_EQ(grids.at("wavelength [m]").NumberOfColumns(), 1);
  ASSERT_EQ(grids.at("wavelength [m]").NumberOfSections(), 156);
  ASSERT_TRUE(grids.at("wavelength [m]").IsConstant());
  ASSERT_NEAR(grids.at("wavelength [m]").edges_(0,0), 120.0*1.0e-9, 1.0e-16);
  ASSERT_NEAR(grids.at("wavelength [m]").edges_(77,0), 311.5*1.0e-9, 1.0e-16);
  ASSERT_NEAR(grids.at("wavelength [m]").edges_(156,0), 735.0*1.0e-9, 1.0e-16);
  ASSERT_NEAR(grids.at("wavelength [m]").mid_points_(0,0), 120.69*1.0e-9, 1.0e-16);
  ASSERT_NEAR(grids.at("wavelength [m]").mid_points_(77,0), 311.9*1.0e-9, 1.0e-16);
  ASSERT_NEAR(grids.at("wavelength [m]").mid_points_(155,0), 729.9*1.0e-9, 1.0e-16);

  // optical depths are unitless
  ASSERT_EQ(accumulated_radiator_states.optical_depth_.Size1(), 156);
  ASSERT_EQ(accumulated_radiator_states.optical_depth_.Size2(), 120);
  ASSERT_EQ(accumulated_radiator_states.optical_depth_.Size3(), 2);
  ASSERT_NEAR(accumulated_radiator_states.optical_depth_(0,0,0), 2896168.4, 1.0e-8);
  ASSERT_NEAR(accumulated_radiator_states.optical_depth_(83,41,0), 4.83152e-4, 1.0e-10);
  ASSERT_NEAR(accumulated_radiator_states.optical_depth_(155,119,0), 6.64755e-10, 1.0e-18);
  ASSERT_NEAR(accumulated_radiator_states.optical_depth_(0,0,1), 2896168.4, 1.0e-3);
  ASSERT_NEAR(accumulated_radiator_states.optical_depth_(83,41,1), 4.83152e-4, 1.0e-8);
  ASSERT_NEAR(accumulated_radiator_states.optical_depth_(155,119,1), 6.64755e-10, 1.0e-18);

  // single scattering albedos are unitless
  ASSERT_EQ(accumulated_radiator_states.single_scattering_albedo_.Size1(), 156);
  ASSERT_EQ(accumulated_radiator_states.single_scattering_albedo_.Size2(), 120);
  ASSERT_EQ(accumulated_radiator_states.single_scattering_albedo_.Size3(), 2);
  ASSERT_NEAR(accumulated_radiator_states.single_scattering_albedo_(0,0,0), 4.78487e-6, 1.0e-12);
  ASSERT_NEAR(accumulated_radiator_states.single_scattering_albedo_(83,41,0), 0.5383450, 1.0e-8);
  ASSERT_NEAR(accumulated_radiator_states.single_scattering_albedo_(155,119,0), 0.997822, 1.0e-8);
  ASSERT_NEAR(accumulated_radiator_states.single_scattering_albedo_(0,0,1), 4.78487e-6, 1.0e-12);
  ASSERT_NEAR(accumulated_radiator_states.single_scattering_albedo_(83,41,1), 0.5383450, 1.0e-8);
  ASSERT_NEAR(accumulated_radiator_states.single_scattering_albedo_(155,119,1), 0.997822, 1.0e-8);

  // asymmetry parameters are unitless
  ASSERT_EQ(accumulated_radiator_states.asymmetry_parameter_.Size1(), 156);
  ASSERT_EQ(accumulated_radiator_states.asymmetry_parameter_.Size2(), 120);
  ASSERT_EQ(accumulated_radiator_states.asymmetry_parameter_.Size3(), 2);
  ASSERT_NEAR(accumulated_radiator_states.asymmetry_parameter_(0,0,0), 2.12438e-2, 1.0e-8);
  ASSERT_NEAR(accumulated_radiator_states.asymmetry_parameter_(83,41,0), 2.13209e-2, 1.0e-8);
  ASSERT_NEAR(accumulated_radiator_states.asymmetry_parameter_(155,119,0), 0.0 , 1.0e-8);
  ASSERT_NEAR(accumulated_radiator_states.asymmetry_parameter_(0,0,1), 2.12438e-2, 1.0e-8);
  ASSERT_NEAR(accumulated_radiator_states.asymmetry_parameter_(83,41,1), 2.13209e-2, 1.0e-8);
  ASSERT_NEAR(accumulated_radiator_states.asymmetry_parameter_(155,119,1), 0.0 , 1.0e-8);

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
  for (int i = 0; i < input.n_wavelengths_; i++) {
    for (int j = 0; j < input.n_levels_; j++) {
      for (int k = 0; k < input.n_columns_; k++) {
        accumulated_radiator_states.optical_depth_(i, j, k) = input.optical_depths_[i * input.n_levels_ * input.n_columns_ + j * input.n_columns_ + k];
        accumulated_radiator_states.single_scattering_albedo_(i, j, k) = input.single_scattering_albedos_[i * input.n_levels_ * input.n_columns_ + j * input.n_columns_ + k];
        accumulated_radiator_states.asymmetry_parameter_(i, j, k) = input.asymmetry_parameters_[i * input.n_levels_ * input.n_columns_ + j * input.n_columns_ + k];
      }
    }
  }
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
