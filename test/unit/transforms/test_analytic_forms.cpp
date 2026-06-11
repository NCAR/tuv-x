// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
//
// Unit tests for the general analytic transform forms, exercising extreme and
// boundary inputs (empty parameter sets, region masking, overflow avoidance,
// values at the band centre and region edges).
#include <tuvx/transforms/analytic_forms.hpp>

#include <gtest/gtest.h>

#include <cmath>
#include <cstddef>
#include <vector>

namespace
{
  // Build a single-column, single-layer state with the given wavelength
  // mid-points (in meters).
  tuvx::AtmosphericState<> StateWithWavelengths(const std::vector<double> &wl_m)
  {
    tuvx::AtmosphericState<> state;
    const auto n = wl_m.size();
    state.wavelength_grid_ = tuvx::Grid<tuvx::Array2D<double>>("m", 1, n);
    for (std::size_t i = 0; i < n; ++i)
    {
      state.wavelength_grid_.mid_points_(i, 0) = wl_m[i];
    }
    state.height_grid_ = tuvx::Grid<tuvx::Array2D<double>>("m", 1, 1);
    state.temperature_ = tuvx::Array2D<double>(1, 1);
    state.pressure_ = tuvx::Array2D<double>(1, 1);
    state.air_density_ = tuvx::Array2D<double>(1, 1);
    return state;
  }
}  // namespace

// ---------------------------------------------------------------------------
// log_normal_bands

TEST(LogNormalBands, EmptyBandsAreZero)
{
  auto state = StateWithWavelengths({ 300e-9, 400e-9 });
  tuvx::Array3D<double> w(2, 1, 1);
  tuvx::log_normal_bands<>({})(state, w);
  EXPECT_DOUBLE_EQ(w(0, 0, 0), 0.0);
  EXPECT_DOUBLE_EQ(w(1, 0, 0), 0.0);
}

TEST(LogNormalBands, PeaksAtCenter)
{
  // ln(center/lambda) = 0 at lambda == center, so the weight equals amplitude.
  auto state = StateWithWavelengths({ 300e-9 });
  tuvx::Array3D<double> w(1, 1, 1);
  tuvx::log_normal_bands<>({ .bands_ = { tuvx::LogNormalBand<>{ .amplitude_ = 5.0, .width_ = 100.0, .center_ = 300e-9 } } })(
      state, w);
  EXPECT_DOUBLE_EQ(w(0, 0, 0), 5.0);
}

TEST(LogNormalBands, SumsMultipleBands)
{
  // Two bands both centred at lambda -> weight is the sum of amplitudes.
  auto state = StateWithWavelengths({ 300e-9 });
  tuvx::Array3D<double> w(1, 1, 1);
  tuvx::log_normal_bands<>({ .bands_ = { tuvx::LogNormalBand<>{ .amplitude_ = 2.0, .width_ = 50.0, .center_ = 300e-9 },
                                         tuvx::LogNormalBand<>{ .amplitude_ = 3.0, .width_ = 50.0, .center_ = 300e-9 } } })(
      state, w);
  EXPECT_DOUBLE_EQ(w(0, 0, 0), 5.0);
}

TEST(LogNormalBands, ZeroOutsideRegion)
{
  auto state = StateWithWavelengths({ 200e-9, 300e-9, 600e-9 });
  tuvx::Array3D<double> w(3, 1, 1);
  tuvx::log_normal_bands<>({ .bands_ = { tuvx::LogNormalBand<>{ .amplitude_ = 1.0, .width_ = 1.0, .center_ = 300e-9 } },
                             .wl_min_ = 250e-9,
                             .wl_max_ = 550e-9 })(state, w);
  EXPECT_DOUBLE_EQ(w(0, 0, 0), 0.0);  // 200 nm: below band
  EXPECT_GT(w(1, 0, 0), 0.0);         // 300 nm: inside
  EXPECT_DOUBLE_EQ(w(2, 0, 0), 0.0);  // 600 nm: above band
}

TEST(LogNormalBands, RegionBoundariesInclusive)
{
  auto state = StateWithWavelengths({ 250e-9, 550e-9 });
  tuvx::Array3D<double> w(2, 1, 1);
  tuvx::log_normal_bands<>({ .bands_ = { tuvx::LogNormalBand<>{ .amplitude_ = 1.0, .width_ = 1.0, .center_ = 400e-9 } },
                             .wl_min_ = 250e-9,
                             .wl_max_ = 550e-9 })(state, w);
  EXPECT_GT(w(0, 0, 0), 0.0);  // exactly at wl_min
  EXPECT_GT(w(1, 0, 0), 0.0);  // exactly at wl_max
}

TEST(LogNormalBands, ExtremeWidthStaysFiniteAndDecays)
{
  // A huge width makes the band extremely narrow; away from centre the
  // exp(-large) underflows to 0 without producing NaN/inf.
  auto state = StateWithWavelengths({ 300e-9, 301e-9 });
  tuvx::Array3D<double> w(2, 1, 1);
  tuvx::log_normal_bands<>(
      { .bands_ = { tuvx::LogNormalBand<>{ .amplitude_ = 1.0, .width_ = 1.0e12, .center_ = 300e-9 } } })(state, w);
  EXPECT_DOUBLE_EQ(w(0, 0, 0), 1.0);  // at centre
  EXPECT_TRUE(std::isfinite(w(1, 0, 0)));
  EXPECT_GE(w(1, 0, 0), 0.0);
  EXPECT_LT(w(1, 0, 0), 1.0e-6);  // far off-peak: underflowed toward 0
}

// ---------------------------------------------------------------------------
// exp_polynomial

TEST(ExpPolynomial, EmptyCoeffsIsExpZero)
{
  // No coefficients -> polynomial is 0 -> exp(0) = 1, times output_scale.
  auto state = StateWithWavelengths({ 300e-9 });
  tuvx::Array3D<double> w(1, 1, 1);
  tuvx::exp_polynomial<>({ .output_scale_ = 7.0 })(state, w);
  EXPECT_DOUBLE_EQ(w(0, 0, 0), 7.0);
}

TEST(ExpPolynomial, ConstantPolynomial)
{
  auto state = StateWithWavelengths({ 300e-9, 400e-9 });
  tuvx::Array3D<double> w(2, 1, 1);
  // exp(2) everywhere, output_scale 1
  tuvx::exp_polynomial<>({ .coefficients_ = { 2.0 } })(state, w);
  EXPECT_DOUBLE_EQ(w(0, 0, 0), std::exp(2.0));
  EXPECT_DOUBLE_EQ(w(1, 0, 0), std::exp(2.0));
}

TEST(ExpPolynomial, LinearWithWavelengthScale)
{
  // poly(x) = 1 + 2*x, evaluated in nm (scale 1e9). At 300 nm: 1 + 2*300 = 601.
  auto state = StateWithWavelengths({ 300e-9 });
  tuvx::Array3D<double> w(1, 1, 1);
  tuvx::exp_polynomial<>({ .coefficients_ = { 1.0, 2.0 }, .wavelength_scale_ = 1.0e9 })(state, w);
  EXPECT_DOUBLE_EQ(w(0, 0, 0), std::exp(1.0 + (2.0 * 300.0)));
}

TEST(ExpPolynomial, OutputScaleApplied)
{
  auto state = StateWithWavelengths({ 300e-9 });
  tuvx::Array3D<double> w(1, 1, 1);
  tuvx::exp_polynomial<>({ .coefficients_ = { 0.0 }, .output_scale_ = 1.0e-4 })(state, w);
  EXPECT_DOUBLE_EQ(w(0, 0, 0), 1.0e-4);  // 1e-4 * exp(0)
}

TEST(ExpPolynomial, HornerMatchesQuadratic)
{
  // poly(x) = c0 + c1*x + c2*x^2 in nm.
  const double c0 = -115.5;
  const double c1 = 0.5307;
  const double c2 = -0.993e-3;
  auto state = StateWithWavelengths({ 300e-9 });
  tuvx::Array3D<double> w(1, 1, 1);
  tuvx::exp_polynomial<>({ .coefficients_ = { c0, c1, c2 },
                           .wavelength_scale_ = 1.0e9,
                           .output_scale_ = 1.0e-4,
                           .wl_min_ = 270e-9,
                           .wl_max_ = 330e-9 })(state, w);
  const double x = 300.0;  // 300 nm
  const double expected = 1.0e-4 * std::exp(c0 + (c1 * x) + (c2 * x * x));
  // Relative tolerance: the form evaluates x = lambda_m * 1e9 by Horner, which
  // differs from this reference expression only by floating-point reassociation.
  EXPECT_NEAR(w(0, 0, 0), expected, std::abs(expected) * 1.0e-12);
}

// The key extreme case: a polynomial that overflows exp() OUTSIDE the active
// region must still yield exactly 0 there -- never inf or NaN -- because the
// formula is not evaluated outside the window.
TEST(ExpPolynomial, NoOverflowOutsideRegion)
{
  // c1 large and positive: exp(huge) would overflow at large wavelengths.
  auto state = StateWithWavelengths({ 100e-9, 300e-9, 900e-9 });
  tuvx::Array3D<double> w(3, 1, 1);
  tuvx::exp_polynomial<>(
      { .coefficients_ = { 0.0, 1.0 }, .wavelength_scale_ = 1.0e9, .wl_min_ = 250e-9, .wl_max_ = 350e-9 })(state, w);
  EXPECT_DOUBLE_EQ(w(0, 0, 0), 0.0);  // 100 nm: outside, exactly zero (not evaluated)
  EXPECT_TRUE(std::isfinite(w(1, 0, 0)));
  EXPECT_DOUBLE_EQ(w(1, 0, 0), std::exp(300.0));  // 300 nm: inside
  EXPECT_DOUBLE_EQ(w(2, 0, 0), 0.0);              // 900 nm: outside; would have overflowed if evaluated
}

TEST(ExpPolynomial, BroadcastsOverHeightAndColumn)
{
  // Two layers, two columns: every (z, col) gets the same per-wavelength value.
  tuvx::AtmosphericState<> state;
  state.wavelength_grid_ = tuvx::Grid<tuvx::Array2D<double>>("m", 1, 1);
  state.wavelength_grid_.mid_points_(0, 0) = 300e-9;
  state.height_grid_ = tuvx::Grid<tuvx::Array2D<double>>("m", 2, 2);
  state.temperature_ = tuvx::Array2D<double>(2, 2);
  state.pressure_ = tuvx::Array2D<double>(2, 2);
  state.air_density_ = tuvx::Array2D<double>(2, 2);

  tuvx::Array3D<double> w(1, 2, 2);
  tuvx::exp_polynomial<>({ .coefficients_ = { 0.5 } })(state, w);
  for (std::size_t z = 0; z < 2; ++z)
  {
    for (std::size_t col = 0; col < 2; ++col)
    {
      EXPECT_DOUBLE_EQ(w(0, z, col), std::exp(0.5));
    }
  }
}
