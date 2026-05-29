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

namespace
{

// Helper to build a simple single-column, single-wavelength test setup.
struct DeltaEddingtonFixture
{
  static constexpr std::size_t N_COLUMNS    = 1;
  static constexpr std::size_t N_LAYERS     = 3;
  static constexpr std::size_t N_WAVELENGTHS = 1;

  tuvx::Grid<tuvx::Array2D<double>> altitude_grid_;
  tuvx::Grid<tuvx::Array2D<double>> wavelength_grid_;
  tuvx::Profile<tuvx::Array2D<double>> surface_albedo_;
  tuvx::ConstituentState<tuvx::Array3D<double>> constituent_state_;
  tuvx::RadiationField<tuvx::RadiationFieldComponents<tuvx::Array3D<double>>> radiation_field_;

  std::map<std::string, tuvx::Grid<tuvx::Array2D<double>>>    grids_;
  std::map<std::string, tuvx::Profile<tuvx::Array2D<double>>> profiles_;

  explicit DeltaEddingtonFixture(double albedo = 0.0)
      : altitude_grid_("m", N_COLUMNS, N_LAYERS),
        wavelength_grid_("m", N_WAVELENGTHS),
        surface_albedo_("1", N_COLUMNS, altitude_grid_),  // dummy grid, overridden below
        constituent_state_(N_COLUMNS, altitude_grid_, wavelength_grid_),
        radiation_field_(N_COLUMNS, altitude_grid_, wavelength_grid_)
  {
    // Altitude edges: ground=0, 1 km, 2 km, 3 km (bottom-to-top)
    for (std::size_t lev = 0; lev <= N_LAYERS; ++lev)
    {
      altitude_grid_.edges_(lev, 0) = static_cast<double>(lev) * 1000.0;
    }
    // Altitude mid-points
    for (std::size_t lay = 0; lay < N_LAYERS; ++lay)
    {
      altitude_grid_.mid_points_(lay, 0) = (static_cast<double>(lay) + 0.5) * 1000.0;
    }

    // Single wavelength band: 300–301 nm (values irrelevant for these tests)
    wavelength_grid_.edges_(0, 0) = 300e-9;
    wavelength_grid_.edges_(1, 0) = 301e-9;
    wavelength_grid_.mid_points_(0, 0) = 300.5e-9;

    // Resize surface albedo to [n_wavelengths × n_columns]
    surface_albedo_ = tuvx::Profile<tuvx::Array2D<double>>(
        "1", N_COLUMNS, wavelength_grid_);
    surface_albedo_.mid_point_values_(0, 0) = albedo;

    grids_["altitude [m]"]    = altitude_grid_;
    grids_["wavelength [m]"]  = wavelength_grid_;
    profiles_["surface albedo [1]"] = surface_albedo_;
  }

  void SetTransparentAtmosphere()
  {
    // Zero optical depth, zero SSA, zero asymmetry
    for (std::size_t wl = 0; wl < N_WAVELENGTHS; ++wl)
    {
      for (std::size_t lay = 0; lay < N_LAYERS; ++lay)
      {
        constituent_state_.optical_depth_(wl, lay, 0)              = 0.0;
        constituent_state_.single_scattering_albedo_(wl, lay, 0)   = 0.0;
        constituent_state_.asymmetry_parameter_(wl, lay, 0)        = 0.0;
      }
    }
  }
};

} // namespace

// For a transparent atmosphere (zero optical depth) with overhead sun (SZA=0)
// and zero surface albedo, the direct irradiance at every level must equal
// mu = cos(0) = 1, and all diffuse components must be zero.
TEST(DeltaEddington, TransparentAtmosphereOverheadSun)
{
  DeltaEddingtonFixture fix(0.0);
  fix.SetTransparentAtmosphere();

  tuvx::DeltaEddington solver;
  const std::vector<double> SZA = { 0.0 };  // overhead sun
  solver.Solve(SZA, fix.grids_, fix.profiles_, fix.constituent_state_, fix.radiation_field_);

  constexpr double TOLERANCE = 1.0e-8;
  constexpr std::size_t N_LEVELS = DeltaEddingtonFixture::N_LAYERS + 1;
  for (std::size_t lev = 0; lev < N_LEVELS; ++lev)
  {
    // Direct irradiance = mu = 1 at all levels (no absorption)
    EXPECT_NEAR(
        fix.radiation_field_.spectral_irradiance_.direct_(0, lev, 0), 1.0, TOLERANCE)
        << "direct irradiance at level " << lev;
    // No diffuse
    EXPECT_NEAR(
        fix.radiation_field_.spectral_irradiance_.upwelling_(0, lev, 0), 0.0, TOLERANCE)
        << "upwelling irradiance at level " << lev;
    EXPECT_NEAR(
        fix.radiation_field_.spectral_irradiance_.downwelling_(0, lev, 0), 0.0, TOLERANCE)
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
  const std::vector<double> SZA = { 0.0 };
  solver.Solve(SZA, fix.grids_, fix.profiles_, fix.constituent_state_, fix.radiation_field_);

  constexpr double TOLERANCE = 1.0e-8;
  // Ground level (index 0 in bottom-to-top convention): direct = mu = 1
  EXPECT_NEAR(
      fix.radiation_field_.spectral_irradiance_.direct_(0, 0, 0), 1.0, TOLERANCE);
}

// Radiation field output dimensions must match what the grids imply.
TEST(DeltaEddington, OutputDimensions)
{
  DeltaEddingtonFixture fix;
  fix.SetTransparentAtmosphere();

  tuvx::DeltaEddington solver;
  const std::vector<double> SZA = { 0.5 };
  solver.Solve(SZA, fix.grids_, fix.profiles_, fix.constituent_state_, fix.radiation_field_);

  // [wavelengths × levels × columns]
  EXPECT_EQ(fix.radiation_field_.spectral_irradiance_.direct_.Size1(),
            DeltaEddingtonFixture::N_WAVELENGTHS);
  EXPECT_EQ(fix.radiation_field_.spectral_irradiance_.direct_.Size2(),
            DeltaEddingtonFixture::N_LAYERS + 1);
  EXPECT_EQ(fix.radiation_field_.spectral_irradiance_.direct_.Size3(),
            DeltaEddingtonFixture::N_COLUMNS);
}

// For a purely absorbing atmosphere (SSA=0, g=0, tau>0) with overhead sun,
// the direct beam must decrease monotonically from TOA to ground.
TEST(DeltaEddington, DirectBeamDecaysWithOpticalDepth)
{
  DeltaEddingtonFixture fix(0.0);

  // Set moderate absorption, no scattering
  for (std::size_t lay = 0; lay < DeltaEddingtonFixture::N_LAYERS; ++lay)
  {
    fix.constituent_state_.optical_depth_(0, lay, 0)            = 0.5;
    fix.constituent_state_.single_scattering_albedo_(0, lay, 0) = 0.0;
    fix.constituent_state_.asymmetry_parameter_(0, lay, 0)      = 0.0;
  }

  tuvx::DeltaEddington solver;
  const std::vector<double> SZA = { 0.0 };
  solver.Solve(SZA, fix.grids_, fix.profiles_, fix.constituent_state_, fix.radiation_field_);

  // TOA level (index n_layers in bottom-to-top) should have the most direct beam
  const double TOA_DIRECT = fix.radiation_field_.spectral_irradiance_.direct_(0, DeltaEddingtonFixture::N_LAYERS, 0);
  EXPECT_NEAR(TOA_DIRECT, 1.0, 1.0e-8);  // mu * exp(0) = 1

  // Ground level (index 0) should have the least
  const double GROUND_DIRECT = fix.radiation_field_.spectral_irradiance_.direct_(0, 0, 0);
  EXPECT_LT(GROUND_DIRECT, TOA_DIRECT);

  // Should decrease monotonically from top to bottom
  for (std::size_t lev = 1; lev <= DeltaEddingtonFixture::N_LAYERS; ++lev)
  {
    const double ABOVE = fix.radiation_field_.spectral_irradiance_.direct_(0, lev, 0);
    const double BELOW = fix.radiation_field_.spectral_irradiance_.direct_(0, lev - 1, 0);
    EXPECT_GE(ABOVE, BELOW) << "direct beam should not increase going down (lev " << lev << ")";
  }
}
