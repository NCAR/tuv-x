// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
//
// Regression tests for DeltaEddington::Solve() against reference CSV data.
//
// Reference CSVs were generated from the Fortran source on the main branch.
// See README.md for the source commit, test atmosphere, and regeneration steps.
//
// Test atmosphere: 3-layer, tau=0.5/layer, ssa=0.9/layer, g=0.85/layer,
//                 surface albedo=0.1, altitudes 0/1/2/3 km.
#include <tuvx/radiative_transfer/constituent.hpp>
#include <tuvx/radiative_transfer/radiation_field.hpp>
#include <tuvx/radiative_transfer/solvers/delta_eddington.hpp>
#include <tuvx/grid.hpp>
#include <tuvx/profile.hpp>

#include <gtest/gtest.h>

#include <cstdlib>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

// ── CSV reader ────────────────────────────────────────────────────────────────

struct ReferenceLevel
{
  std::size_t level;
  double      direct_irradiance;
  double      upwelling_irradiance;
  double      downwelling_irradiance;
  double      direct_actinic_flux;
  double      upwelling_actinic_flux;
  double      downwelling_actinic_flux;
};

static std::vector<ReferenceLevel> LoadReference(const std::string& path)
{
  std::ifstream file(path);
  EXPECT_TRUE(file.is_open()) << "Cannot open reference file: " << path;

  std::vector<ReferenceLevel> rows;
  std::string line;
  while (std::getline(file, line))
  {
    if (line.empty() || line[0] == '#' || line[0] == '"')
    {
      continue;
    }
    if (line.find("level") != std::string::npos)
    {
      continue;  // header row
    }
    std::istringstream ss(line);
    std::string token;
    ReferenceLevel row{};
    std::getline(ss, token, ','); row.level                    = std::stoul(token);
    std::getline(ss, token, ','); row.direct_irradiance        = std::stod(token);
    std::getline(ss, token, ','); row.upwelling_irradiance     = std::stod(token);
    std::getline(ss, token, ','); row.downwelling_irradiance   = std::stod(token);
    std::getline(ss, token, ','); row.direct_actinic_flux      = std::stod(token);
    std::getline(ss, token, ','); row.upwelling_actinic_flux   = std::stod(token);
    std::getline(ss, token, ','); row.downwelling_actinic_flux = std::stod(token);
    rows.push_back(row);
  }
  return rows;
}

// ── Test fixture ──────────────────────────────────────────────────────────────

struct RegressionAtmosphere
{
  static constexpr std::size_t n_columns    = 1;
  static constexpr std::size_t n_layers     = 3;
  static constexpr std::size_t n_wavelengths = 1;

  // Optical properties: top-to-bottom
  static constexpr double optical_depth_per_layer       = 0.5;
  static constexpr double single_scattering_albedo_value = 0.9;
  static constexpr double asymmetry_parameter_value      = 0.85;
  static constexpr double surface_albedo_value           = 0.1;

  tuvx::Grid<tuvx::Array2D<double>>    altitude_grid;
  tuvx::Grid<tuvx::Array2D<double>>    wavelength_grid;
  tuvx::Profile<tuvx::Array2D<double>> surface_albedo;
  tuvx::ConstituentState<tuvx::Array3D<double>> constituent_state;
  tuvx::RadiationField<tuvx::RadiationFieldComponents<tuvx::Array3D<double>>> radiation_field;

  std::map<std::string, tuvx::Grid<tuvx::Array2D<double>>>    grids;
  std::map<std::string, tuvx::Profile<tuvx::Array2D<double>>> profiles;

  RegressionAtmosphere()
      : altitude_grid("m", n_columns, n_layers),
        wavelength_grid("m", n_wavelengths),
        surface_albedo("1", n_columns, altitude_grid),
        constituent_state(n_columns, altitude_grid, wavelength_grid),
        radiation_field(n_columns, altitude_grid, wavelength_grid)
  {
    // Altitude edges: 0, 1, 2, 3 km (bottom-to-top)
    for (std::size_t lev = 0; lev <= n_layers; ++lev)
    {
      altitude_grid.edges_(lev, 0) = static_cast<double>(lev) * 1000.0;
    }
    for (std::size_t lay = 0; lay < n_layers; ++lay)
    {
      altitude_grid.mid_points_(lay, 0) = (static_cast<double>(lay) + 0.5) * 1000.0;
    }

    wavelength_grid.edges_(0, 0)     = 300e-9;
    wavelength_grid.edges_(1, 0)     = 301e-9;
    wavelength_grid.mid_points_(0, 0) = 300.5e-9;

    surface_albedo = tuvx::Profile<tuvx::Array2D<double>>("1", n_columns, wavelength_grid);
    surface_albedo.mid_point_values_(0, 0) = surface_albedo_value;

    // Uniform optical properties, top-to-bottom
    for (std::size_t lay = 0; lay < n_layers; ++lay)
    {
      constituent_state.optical_depth_(0, lay, 0)            = optical_depth_per_layer;
      constituent_state.single_scattering_albedo_(0, lay, 0) = single_scattering_albedo_value;
      constituent_state.asymmetry_parameter_(0, lay, 0)      = asymmetry_parameter_value;
    }

    grids["altitude [m]"]       = altitude_grid;
    grids["wavelength [m]"]     = wavelength_grid;
    profiles["surface albedo [1]"] = surface_albedo;
  }
};

// ── Helper ────────────────────────────────────────────────────────────────────

static std::string ReferenceFile(double sza_rad)
{
  // Locate the reference directory relative to this source file.
  // CMAKE_CURRENT_SOURCE_DIR is injected via a compile definition.
  const char* src_dir = std::getenv("TUVX_REGRESSION_REFERENCE_DIR");
  std::string dir = src_dir ? src_dir : TUVX_REGRESSION_REFERENCE_DIR;
  char buf[64];
  std::snprintf(buf, sizeof(buf), "delta_eddington_sza_%.4f.csv", sza_rad);
  return dir + "/" + buf;
}

// ── Tests ─────────────────────────────────────────────────────────────────────

static void RunRegressionTest(double sza_rad, double tolerance)
{
  const std::string path = ReferenceFile(sza_rad);
  const auto reference   = LoadReference(path);
  ASSERT_FALSE(reference.empty()) << "No rows loaded from " << path;

  RegressionAtmosphere atm;
  tuvx::DeltaEddington solver;
  const std::vector<double> sza_vec = { sza_rad };
  solver.Solve(sza_vec, atm.grids, atm.profiles, atm.constituent_state, atm.radiation_field);

  for (const auto& ref : reference)
  {
    const std::size_t lev = ref.level;
    EXPECT_NEAR(atm.radiation_field.spectral_irradiance_.direct_(0, lev, 0),
                ref.direct_irradiance, tolerance)
        << "direct_irradiance mismatch at level " << lev << " sza=" << sza_rad;
    EXPECT_NEAR(atm.radiation_field.spectral_irradiance_.upwelling_(0, lev, 0),
                ref.upwelling_irradiance, tolerance)
        << "upwelling_irradiance mismatch at level " << lev << " sza=" << sza_rad;
    EXPECT_NEAR(atm.radiation_field.spectral_irradiance_.downwelling_(0, lev, 0),
                ref.downwelling_irradiance, tolerance)
        << "downwelling_irradiance mismatch at level " << lev << " sza=" << sza_rad;
    EXPECT_NEAR(atm.radiation_field.actinic_flux_.direct_(0, lev, 0),
                ref.direct_actinic_flux, tolerance)
        << "direct_actinic_flux mismatch at level " << lev << " sza=" << sza_rad;
    EXPECT_NEAR(atm.radiation_field.actinic_flux_.upwelling_(0, lev, 0),
                ref.upwelling_actinic_flux, tolerance)
        << "upwelling_actinic_flux mismatch at level " << lev << " sza=" << sza_rad;
    EXPECT_NEAR(atm.radiation_field.actinic_flux_.downwelling_(0, lev, 0),
                ref.downwelling_actinic_flux, tolerance)
        << "downwelling_actinic_flux mismatch at level " << lev << " sza=" << sza_rad;
  }
}

TEST(DeltaEddingtonRegression, Sza0_0)
{
  RunRegressionTest(0.0, 1.0e-10);
}

TEST(DeltaEddingtonRegression, Sza0_5)
{
  RunRegressionTest(0.5, 1.0e-10);
}

TEST(DeltaEddingtonRegression, Sza1_0)
{
  RunRegressionTest(1.0, 1.0e-10);
}
