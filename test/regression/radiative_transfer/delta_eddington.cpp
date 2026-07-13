// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
//
// Regression tests for DeltaEddington::Solve() against reference CSV data.
//
// Test atmosphere: 3-layer, tau=0.5/layer, ssa=0.9/layer, g=0.85/layer,
//                 surface albedo=0.1, altitudes 0/1/2/3 km.
#include <tuvx/grid.hpp>
#include <tuvx/profile.hpp>
#include <tuvx/radiative_transfer/constituent.hpp>
#include <tuvx/radiative_transfer/radiation_field.hpp>
#include <tuvx/radiative_transfer/solvers/delta_eddington.hpp>

#include <gtest/gtest.h>

#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <map>
#include <sstream>
#include <string>
#include <vector>

// ── CSV reader ────────────────────────────────────────────────────────────────

namespace
{

  struct ReferenceLevel
  {
    std::size_t level_ = 0;
    double direct_irradiance_ = 0.0;
    double upwelling_irradiance_ = 0.0;
    double downwelling_irradiance_ = 0.0;
    double direct_actinic_flux_ = 0.0;
    double upwelling_actinic_flux_ = 0.0;
    double downwelling_actinic_flux_ = 0.0;
  };

  std::vector<ReferenceLevel> LoadReference(const std::string& path)
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
      std::getline(ss, token, ',');
      row.level_ = std::stoul(token);
      std::getline(ss, token, ',');
      row.direct_irradiance_ = std::stod(token);
      std::getline(ss, token, ',');
      row.upwelling_irradiance_ = std::stod(token);
      std::getline(ss, token, ',');
      row.downwelling_irradiance_ = std::stod(token);
      std::getline(ss, token, ',');
      row.direct_actinic_flux_ = std::stod(token);
      std::getline(ss, token, ',');
      row.upwelling_actinic_flux_ = std::stod(token);
      std::getline(ss, token, ',');
      row.downwelling_actinic_flux_ = std::stod(token);
      rows.push_back(row);
    }
    return rows;
  }

  // ── Test fixture ──────────────────────────────────────────────────────────────

  struct RegressionAtmosphere
  {
    static constexpr std::size_t N_COLUMNS = 1;
    static constexpr std::size_t N_LAYERS = 3;
    static constexpr std::size_t N_WAVELENGTHS = 1;

    // Optical properties: top-to-bottom
    static constexpr double OPTICAL_DEPTH_PER_LAYER = 0.5;
    static constexpr double SINGLE_SCATTERING_ALBEDO_VALUE = 0.9;
    static constexpr double ASYMMETRY_PARAMETER_VALUE = 0.85;
    static constexpr double SURFACE_ALBEDO_VALUE = 0.1;

    tuvx::Grid<tuvx::Array2D<double>> altitude_grid_;
    tuvx::Grid<tuvx::Array2D<double>> wavelength_grid_;
    tuvx::Profile<tuvx::Array2D<double>> surface_albedo_;
    tuvx::ConstituentState<tuvx::Array3D<double>> constituent_state_;
    tuvx::RadiationField<tuvx::RadiationFieldComponents<tuvx::Array3D<double>>> radiation_field_;

    std::map<std::string, tuvx::Grid<tuvx::Array2D<double>>> grids_;
    std::map<std::string, tuvx::Profile<tuvx::Array2D<double>>> profiles_;

    RegressionAtmosphere()
        : altitude_grid_("m", N_COLUMNS, N_LAYERS),
          wavelength_grid_("m", N_WAVELENGTHS),
          surface_albedo_("1", N_COLUMNS, altitude_grid_),
          constituent_state_(N_COLUMNS, altitude_grid_, wavelength_grid_),
          radiation_field_(N_COLUMNS, altitude_grid_, wavelength_grid_)
    {
      // Altitude edges: 0, 1, 2, 3 km (bottom-to-top)
      for (std::size_t lev = 0; lev <= N_LAYERS; ++lev)
      {
        altitude_grid_.edges_(lev, 0) = static_cast<double>(lev) * 1000.0;
      }
      for (std::size_t lay = 0; lay < N_LAYERS; ++lay)
      {
        altitude_grid_.mid_points_(lay, 0) = (static_cast<double>(lay) + 0.5) * 1000.0;
      }

      wavelength_grid_.edges_(0, 0) = 300e-9;
      wavelength_grid_.edges_(1, 0) = 301e-9;
      wavelength_grid_.mid_points_(0, 0) = 300.5e-9;

      surface_albedo_ = tuvx::Profile<tuvx::Array2D<double>>("1", N_COLUMNS, wavelength_grid_);
      surface_albedo_.mid_point_values_(0, 0) = SURFACE_ALBEDO_VALUE;

      // Uniform optical properties, top-to-bottom
      for (std::size_t lay = 0; lay < N_LAYERS; ++lay)
      {
        constituent_state_.optical_depth_(0, lay, 0) = OPTICAL_DEPTH_PER_LAYER;
        constituent_state_.single_scattering_albedo_(0, lay, 0) = SINGLE_SCATTERING_ALBEDO_VALUE;
        constituent_state_.asymmetry_parameter_(0, lay, 0) = ASYMMETRY_PARAMETER_VALUE;
      }

      grids_["altitude [m]"] = altitude_grid_;
      grids_["wavelength [m]"] = wavelength_grid_;
      profiles_["surface albedo [1]"] = surface_albedo_;
    }
  };

  // ── Helper ────────────────────────────────────────────────────────────────────

  std::string ReferenceFile(double sza_rad)
  {
    const char* src_dir = std::getenv("TUVX_REGRESSION_REFERENCE_DIR");
    const std::string DIR = (src_dir != nullptr) ? src_dir : TUVX_REGRESSION_REFERENCE_DIR;
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(4) << sza_rad;
    return DIR + "/delta_eddington_sza_" + oss.str() + ".csv";
  }

  void CheckReferenceLevel(
      const tuvx::RadiationField<tuvx::RadiationFieldComponents<tuvx::Array3D<double>>>& field,
      const ReferenceLevel& ref,
      double tolerance,
      double sza_rad)
  {
    const std::size_t LEV = ref.level_;
    EXPECT_NEAR(field.spectral_irradiance_.direct_(0, LEV, 0), ref.direct_irradiance_, tolerance)
        << "direct_irradiance mismatch at level " << LEV << " sza=" << sza_rad;
    EXPECT_NEAR(field.spectral_irradiance_.upwelling_(0, LEV, 0), ref.upwelling_irradiance_, tolerance)
        << "upwelling_irradiance mismatch at level " << LEV << " sza=" << sza_rad;
    EXPECT_NEAR(field.spectral_irradiance_.downwelling_(0, LEV, 0), ref.downwelling_irradiance_, tolerance)
        << "downwelling_irradiance mismatch at level " << LEV << " sza=" << sza_rad;
    EXPECT_NEAR(field.actinic_flux_.direct_(0, LEV, 0), ref.direct_actinic_flux_, tolerance)
        << "direct_actinic_flux mismatch at level " << LEV << " sza=" << sza_rad;
    EXPECT_NEAR(field.actinic_flux_.upwelling_(0, LEV, 0), ref.upwelling_actinic_flux_, tolerance)
        << "upwelling_actinic_flux mismatch at level " << LEV << " sza=" << sza_rad;
    EXPECT_NEAR(field.actinic_flux_.downwelling_(0, LEV, 0), ref.downwelling_actinic_flux_, tolerance)
        << "downwelling_actinic_flux mismatch at level " << LEV << " sza=" << sza_rad;
  }

  // ── Tests ─────────────────────────────────────────────────────────────────────

  void RunRegressionTest(double sza_rad, double tolerance)
  {
    const std::string PATH = ReferenceFile(sza_rad);
    const auto REFERENCE = LoadReference(PATH);
    ASSERT_FALSE(REFERENCE.empty()) << "No rows loaded from " << PATH;

    RegressionAtmosphere atm;
    tuvx::DeltaEddington solver;
    const std::vector<double> SZA_VEC = { sza_rad };
    solver.Solve(SZA_VEC, atm.grids_, atm.profiles_, atm.constituent_state_, atm.radiation_field_);

    for (const auto& ref : REFERENCE)
    {
      CheckReferenceLevel(atm.radiation_field_, ref, tolerance, sza_rad);
    }
  }

}  // namespace

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
