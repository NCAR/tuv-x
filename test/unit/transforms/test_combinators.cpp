// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#include <tuvx/transforms/combinators.hpp>
#include <tuvx/transforms/factories.hpp>

#include <gtest/gtest.h>

namespace
{
  tuvx::AtmosphericState<> make_state()
  {
    tuvx::AtmosphericState<> state;

    state.wavelength_grid_ = tuvx::Grid<tuvx::Array2D<double>>("m", 1, 3);
    state.wavelength_grid_.edges_(0, 0) = 200e-9;
    state.wavelength_grid_.edges_(1, 0) = 250e-9;
    state.wavelength_grid_.edges_(2, 0) = 300e-9;
    state.wavelength_grid_.edges_(3, 0) = 400e-9;
    state.wavelength_grid_.mid_points_(0, 0) = 225e-9;
    state.wavelength_grid_.mid_points_(1, 0) = 275e-9;
    state.wavelength_grid_.mid_points_(2, 0) = 350e-9;

    state.height_grid_ = tuvx::Grid<tuvx::Array2D<double>>("m", 2, 2);

    state.temperature_ = tuvx::Array2D<double>(2, 2);
    state.temperature_(0, 0) = 250.0;
    state.temperature_(1, 0) = 270.0;
    state.temperature_(0, 1) = 260.0;
    state.temperature_(1, 1) = 280.0;

    state.pressure_ = tuvx::Array2D<double>(2, 2);
    state.air_density_ = tuvx::Array2D<double>(2, 2);
    state.air_density_(0, 0) = 1.0e-3;
    state.air_density_(1, 0) = 8.0e-4;
    state.air_density_(0, 1) = 9.0e-4;
    state.air_density_(1, 1) = 7.0e-4;

    return state;
  }

  tuvx::Array3D<double> make_weights()
  {
    return { 3, 2, 2 };
  }
}  // namespace

// ---------------------------------------------------------------------------

TEST(Combinators, Multiply)
{
  auto state = make_state();
  auto weights = make_weights();

  auto tf = tuvx::multiply(tuvx::constant(2.0), tuvx::constant(3.0));
  tf(state, weights);

  for (const auto &w : weights)
  {
    EXPECT_DOUBLE_EQ(w, 6.0);
  }
}

TEST(Combinators, Add)
{
  auto state = make_state();
  auto weights = make_weights();

  auto tf = tuvx::add(tuvx::constant(4.0), tuvx::constant(1.5));
  tf(state, weights);

  for (const auto &w : weights)
  {
    EXPECT_DOUBLE_EQ(w, 5.5);
  }
}

TEST(Combinators, Scale)
{
  auto state = make_state();
  auto weights = make_weights();

  auto tf = tuvx::scale(5.0, tuvx::constant(2.0));
  tf(state, weights);

  for (const auto &w : weights)
  {
    EXPECT_DOUBLE_EQ(w, 10.0);
  }
}

// in_region: wl mid-points 225nm, 275nm, 350nm
// region [250nm, 300nm] -> only bin 1 (275nm) is inside
TEST(Combinators, InRegion)
{
  auto state = make_state();
  auto weights = make_weights();

  auto tf = tuvx::in_region(250e-9, 300e-9, tuvx::constant(1.0));
  tf(state, weights);

  for (std::size_t z = 0; z < 2; ++z)
  {
    for (std::size_t col = 0; col < 2; ++col)
    {
      EXPECT_DOUBLE_EQ(weights(0, z, col), 0.0);
      EXPECT_DOUBLE_EQ(weights(1, z, col), 1.0);
      EXPECT_DOUBLE_EQ(weights(2, z, col), 0.0);
    }
  }
}

TEST(Combinators, ClampHigh)
{
  auto state = make_state();
  auto weights = make_weights();

  auto tf = tuvx::clamp(tuvx::constant(5.0), 0.0, 3.0);
  tf(state, weights);

  for (const auto &w : weights)
  {
    EXPECT_DOUBLE_EQ(w, 3.0);
  }
}

TEST(Combinators, ClampLow)
{
  auto state = make_state();
  auto weights = make_weights();

  auto tf = tuvx::clamp(tuvx::constant(-1.0), 0.0, 3.0);
  tf(state, weights);

  for (const auto &w : weights)
  {
    EXPECT_DOUBLE_EQ(w, 0.0);
  }
}

// piecewise: 3 regions cover the 3 wavelength bins with distinct values
TEST(Combinators, Piecewise)
{
  auto state = make_state();
  auto weights = make_weights();

  std::vector<tuvx::PiecewiseRegion<>> regions;
  regions.push_back({ .wl_min_ = 200e-9, .wl_max_ = 249e-9, .transform_ = tuvx::constant(10.0) });
  regions.push_back({ .wl_min_ = 250e-9, .wl_max_ = 310e-9, .transform_ = tuvx::constant(20.0) });
  regions.push_back({ .wl_min_ = 310e-9, .wl_max_ = 400e-9, .transform_ = tuvx::constant(30.0) });

  auto tf = tuvx::piecewise(std::move(regions));
  tf(state, weights);

  for (std::size_t z = 0; z < 2; ++z)
  {
    for (std::size_t col = 0; col < 2; ++col)
    {
      EXPECT_DOUBLE_EQ(weights(0, z, col), 10.0);
      EXPECT_DOUBLE_EQ(weights(1, z, col), 20.0);
      EXPECT_DOUBLE_EQ(weights(2, z, col), 30.0);
    }
  }
}

TEST(Combinators, PiecewiseGapIsZero)
{
  auto state = make_state();
  auto weights = make_weights();

  // Only covers bins 0 and 2  -  bin 1 (275nm) gets zero
  std::vector<tuvx::PiecewiseRegion<>> regions;
  regions.push_back({ .wl_min_ = 200e-9, .wl_max_ = 249e-9, .transform_ = tuvx::constant(1.0) });
  regions.push_back({ .wl_min_ = 310e-9, .wl_max_ = 400e-9, .transform_ = tuvx::constant(2.0) });

  auto tf = tuvx::piecewise(std::move(regions));
  tf(state, weights);

  for (std::size_t z = 0; z < 2; ++z)
  {
    for (std::size_t col = 0; col < 2; ++col)
    {
      EXPECT_DOUBLE_EQ(weights(0, z, col), 1.0);
      EXPECT_DOUBLE_EQ(weights(1, z, col), 0.0);
      EXPECT_DOUBLE_EQ(weights(2, z, col), 2.0);
    }
  }
}

// override_band: schumann-runge is [175-205nm]; test wavelengths are outside
TEST(Combinators, OverrideBandNoOverlap)
{
  auto state = make_state();
  auto weights = make_weights();

  auto tf = tuvx::override_band("schumann-runge", 99.0, tuvx::constant(1.0));
  tf(state, weights);

  for (const auto &w : weights)
  {
    EXPECT_DOUBLE_EQ(w, 1.0);
  }
}

// override_band: lyman-alpha [121.50nm, 121.65nm]; build state with bin inside band
TEST(Combinators, OverrideBandHit)
{
  tuvx::AtmosphericState<> state;
  state.wavelength_grid_ = tuvx::Grid<tuvx::Array2D<double>>("m", 1, 1);
  state.wavelength_grid_.edges_(0, 0) = 1.21e-7;
  state.wavelength_grid_.edges_(1, 0) = 1.22e-7;
  state.wavelength_grid_.mid_points_(0, 0) = 1.215e-7;

  state.height_grid_ = tuvx::Grid<tuvx::Array2D<double>>("m", 1, 1);
  state.temperature_ = tuvx::Array2D<double>(1, 1);
  state.temperature_(0, 0) = 250.0;
  state.pressure_ = tuvx::Array2D<double>(1, 1);
  state.air_density_ = tuvx::Array2D<double>(1, 1);
  state.air_density_(0, 0) = 1.0e-3;

  tuvx::Array3D<double> weights(1, 1, 1);
  auto tf = tuvx::override_band("lyman-alpha", 42.0, tuvx::constant(1.0));
  tf(state, weights);

  EXPECT_DOUBLE_EQ(weights(0, 0, 0), 42.0);
}

TEST(Combinators, OverrideBandInvalidName)
{
  EXPECT_THROW(tuvx::named_band("nonexistent"), std::invalid_argument);
}

// Composition: add(base, in_region(...)) and multiply(base, in_region(...))
TEST(Combinators, Composition)
{
  auto state = make_state();

  // add: base(1e-22) + in_region(250-300, 1e-22) -> 1e-22 except bin1 gets 2e-22
  auto weights_add = make_weights();
  auto tf_add = tuvx::add(tuvx::constant(1.0e-22), tuvx::in_region(250e-9, 300e-9, tuvx::constant(1.0e-22)));
  tf_add(state, weights_add);

  for (std::size_t z = 0; z < 2; ++z)
  {
    for (std::size_t col = 0; col < 2; ++col)
    {
      EXPECT_DOUBLE_EQ(weights_add(0, z, col), 1.0e-22);
      EXPECT_DOUBLE_EQ(weights_add(1, z, col), 2.0e-22);
      EXPECT_DOUBLE_EQ(weights_add(2, z, col), 1.0e-22);
    }
  }

  // multiply: 3 x in_region(250-300, 2) -> 0 outside, 6 inside
  auto weights_mul = make_weights();
  auto tf_mul = tuvx::multiply(tuvx::constant(3.0), tuvx::in_region(250e-9, 300e-9, tuvx::constant(2.0)));
  tf_mul(state, weights_mul);

  for (std::size_t z = 0; z < 2; ++z)
  {
    for (std::size_t col = 0; col < 2; ++col)
    {
      EXPECT_DOUBLE_EQ(weights_mul(0, z, col), 0.0);
      EXPECT_DOUBLE_EQ(weights_mul(1, z, col), 6.0);
      EXPECT_DOUBLE_EQ(weights_mul(2, z, col), 0.0);
    }
  }
}
