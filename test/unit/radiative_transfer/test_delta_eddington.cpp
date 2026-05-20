// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#include <tuvx/radiative_transfer/constituent.hpp>
#include <tuvx/radiative_transfer/radiation_field.hpp>
#include <tuvx/radiative_transfer/solvers/delta_eddington.hpp>
#include <tuvx/grid.hpp>
#include <tuvx/profile.hpp>

#include <gtest/gtest.h>

#include <cmath>
#include <map>
#include <string>
#include <vector>

// Helper to build a simple single-column, single-wavelength test setup.
// Returns grids, profiles, constituent state, and radiation field.
struct DeltaEddingtonFixture
{
  static constexpr std::size_t kNColumns    = 1;
  static constexpr std::size_t kNLayers     = 3;
  static constexpr std::size_t kNWavelengths = 1;

  tuvx::Grid<tuvx::Array2D<double>> altitude_grid;
  tuvx::Grid<tuvx::Array2D<double>> wavelength_grid;
  tuvx::Profile<tuvx::Array2D<double>> surface_albedo;
  tuvx::ConstituentState<tuvx::Array3D<double>> constituent_state;
  tuvx::RadiationField<tuvx::RadiationFieldComponents<tuvx::Array3D<double>>> radiation_field;

  std::map<std::string, tuvx::Grid<tuvx::Array2D<double>>>    grids;
  std::map<std::string, tuvx::Profile<tuvx::Array2D<double>>> profiles;

  explicit DeltaEddingtonFixture(double albedo = 0.0)
      : altitude_grid("m", kNColumns, kNLayers),
        wavelength_grid("m", kNWavelengths),
        surface_albedo("1", kNColumns, altitude_grid),  // dummy grid, overridden below
        constituent_state(kNColumns, altitude_grid, wavelength_grid),
        radiation_field(kNColumns, altitude_grid, wavelength_grid)
  {
    // Altitude edges: ground=0, 1 km, 2 km, 3 km (bottom-to-top)
    for (std::size_t lev = 0; lev <= kNLayers; ++lev)
    {
      altitude_grid.edges_(lev, 0) = static_cast<double>(lev) * 1000.0;
    }
    // Altitude mid-points
    for (std::size_t lay = 0; lay < kNLayers; ++lay)
    {
      altitude_grid.mid_points_(lay, 0) = (static_cast<double>(lay) + 0.5) * 1000.0;
    }

    // Single wavelength band: 300–301 nm (values irrelevant for these tests)
    wavelength_grid.edges_(0, 0) = 300e-9;
    wavelength_grid.edges_(1, 0) = 301e-9;
    wavelength_grid.mid_points_(0, 0) = 300.5e-9;

    // Resize surface albedo to [n_wavelengths × n_columns]
    surface_albedo = tuvx::Profile<tuvx::Array2D<double>>(
        "1", kNColumns, wavelength_grid);
    surface_albedo.mid_point_values_(0, 0) = albedo;

    grids["altitude [m]"]    = altitude_grid;
    grids["wavelength [m]"]  = wavelength_grid;
    profiles["surface albedo [1]"] = surface_albedo;
  }

  void SetTransparentAtmosphere()
  {
    // Zero optical depth, zero SSA, zero asymmetry
    for (std::size_t wl = 0; wl < kNWavelengths; ++wl)
    {
      for (std::size_t lay = 0; lay < kNLayers; ++lay)
      {
        constituent_state.optical_depth_(wl, lay, 0)              = 0.0;
        constituent_state.single_scattering_albedo_(wl, lay, 0)   = 0.0;
        constituent_state.asymmetry_parameter_(wl, lay, 0)        = 0.0;
      }
    }
  }
};

// For a transparent atmosphere (zero optical depth) with overhead sun (SZA=0)
// and zero surface albedo, the direct irradiance at every level must equal
// mu = cos(0) = 1, and all diffuse components must be zero.
TEST(DeltaEddington, TransparentAtmosphereOverheadSun)
{
  DeltaEddingtonFixture fix(0.0);
  fix.SetTransparentAtmosphere();

  tuvx::DeltaEddington solver;
  const std::vector<double> sza = { 0.0 };  // overhead sun
  solver.Solve(sza, fix.grids, fix.profiles, fix.constituent_state, fix.radiation_field);

  constexpr double kTol = 1.0e-8;
  constexpr std::size_t n_levels = DeltaEddingtonFixture::kNLayers + 1;
  for (std::size_t lev = 0; lev < n_levels; ++lev)
  {
    // Direct irradiance = mu = 1 at all levels (no absorption)
    EXPECT_NEAR(
        fix.radiation_field.spectral_irradiance_.direct_(0, lev, 0), 1.0, kTol)
        << "direct irradiance at level " << lev;
    // No diffuse
    EXPECT_NEAR(
        fix.radiation_field.spectral_irradiance_.upwelling_(0, lev, 0), 0.0, kTol)
        << "upwelling irradiance at level " << lev;
    EXPECT_NEAR(
        fix.radiation_field.spectral_irradiance_.downwelling_(0, lev, 0), 0.0, kTol)
        << "downwelling irradiance at level " << lev;
  }
}

// For a transparent atmosphere with overhead sun and non-zero surface albedo,
// the direct beam still reaches the ground unattenuated.
TEST(DeltaEddington, TransparentAtmosphereWithAlbedo)
{
  DeltaEddingtonFixture fix(0.1);
  fix.SetTransparentAtmosphere();

  tuvx::DeltaEddington solver;
  const std::vector<double> sza = { 0.0 };
  solver.Solve(sza, fix.grids, fix.profiles, fix.constituent_state, fix.radiation_field);

  constexpr double kTol = 1.0e-8;
  // Ground level (index 0 in bottom-to-top convention): direct = mu = 1
  EXPECT_NEAR(
      fix.radiation_field.spectral_irradiance_.direct_(0, 0, 0), 1.0, kTol);
}

// Radiation field output dimensions must match what the grids imply.
TEST(DeltaEddington, OutputDimensions)
{
  DeltaEddingtonFixture fix;
  fix.SetTransparentAtmosphere();

  tuvx::DeltaEddington solver;
  const std::vector<double> sza = { 0.5 };
  solver.Solve(sza, fix.grids, fix.profiles, fix.constituent_state, fix.radiation_field);

  // [wavelengths × levels × columns]
  EXPECT_EQ(fix.radiation_field.spectral_irradiance_.direct_.Size1(),
            DeltaEddingtonFixture::kNWavelengths);
  EXPECT_EQ(fix.radiation_field.spectral_irradiance_.direct_.Size2(),
            DeltaEddingtonFixture::kNLayers + 1);
  EXPECT_EQ(fix.radiation_field.spectral_irradiance_.direct_.Size3(),
            DeltaEddingtonFixture::kNColumns);
}

// For a purely absorbing atmosphere (SSA=0, g=0, tau>0) with overhead sun,
// the direct beam must decrease monotonically from TOA to ground.
TEST(DeltaEddington, DirectBeamDecaysWithOpticalDepth)
{
  DeltaEddingtonFixture fix(0.0);

  // Set moderate absorption, no scattering
  for (std::size_t lay = 0; lay < DeltaEddingtonFixture::kNLayers; ++lay)
  {
    fix.constituent_state.optical_depth_(0, lay, 0)            = 0.5;
    fix.constituent_state.single_scattering_albedo_(0, lay, 0) = 0.0;
    fix.constituent_state.asymmetry_parameter_(0, lay, 0)      = 0.0;
  }

  tuvx::DeltaEddington solver;
  const std::vector<double> sza = { 0.0 };
  solver.Solve(sza, fix.grids, fix.profiles, fix.constituent_state, fix.radiation_field);

  // TOA level (index n_layers in bottom-to-top) should have the most direct beam
  constexpr std::size_t n_layers = DeltaEddingtonFixture::kNLayers;
  const double toa_direct = fix.radiation_field.spectral_irradiance_.direct_(0, n_layers, 0);
  EXPECT_NEAR(toa_direct, 1.0, 1.0e-8);  // mu * exp(0) = 1

  // Ground level (index 0) should have the least
  const double ground_direct = fix.radiation_field.spectral_irradiance_.direct_(0, 0, 0);
  EXPECT_LT(ground_direct, toa_direct);

  // Should decrease monotonically from top to bottom
  for (std::size_t lev = 1; lev <= n_layers; ++lev)
  {
    const double above = fix.radiation_field.spectral_irradiance_.direct_(0, lev, 0);
    const double below = fix.radiation_field.spectral_irradiance_.direct_(0, lev - 1, 0);
    EXPECT_GE(above, below) << "direct beam should not increase going down (lev " << lev << ")";
  }
}
