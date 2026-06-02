// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#include <tuvx/transforms/factories.hpp>

#include <gtest/gtest.h>

#include <cmath>

namespace
{
  // Test atmosphere: 3 wavelength bins, 2 height layers, 2 columns.
  //   wl mid-points: 225nm, 275nm, 350nm  (edges: 200, 250, 300, 400 nm)
  //   temperature:   col0: [250, 270],  col1: [260, 280]  (K)
  //   air density:   col0: [1.0e-3, 8.0e-4],  col1: [9.0e-4, 7.0e-4]  (mol/m3)
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
    state.pressure_(0, 0) = 1.0e5;
    state.pressure_(1, 0) = 0.5e5;
    state.pressure_(0, 1) = 1.0e5;
    state.pressure_(1, 1) = 0.5e5;

    state.air_density_ = tuvx::Array2D<double>(2, 2);
    state.air_density_(0, 0) = 1.0e-3;
    state.air_density_(1, 0) = 8.0e-4;
    state.air_density_(0, 1) = 9.0e-4;
    state.air_density_(1, 1) = 7.0e-4;

    return state;
  }

  // weights shape: [n_wl=3, n_z=2, n_col=2]
  tuvx::Array3D<double> make_weights()
  {
    return { 3, 2, 2 };
  }
}  // namespace

// ---------------------------------------------------------------------------

TEST(Factories, Constant)
{
  auto state = make_state();
  auto weights = make_weights();

  const auto tf = tuvx::constant(3.14);
  tf(state, weights);

  for (const auto &w : weights)
  {
    EXPECT_DOUBLE_EQ(w, 3.14);
  }
}

TEST(Factories, WrapAnalyticWavelengthOnly)
{
  auto state = make_state();
  auto weights = make_weights();

  // f(lambda) = lambda in nm -> 225, 275, 350
  auto tf = tuvx::wrap_analytic([](double lambda) { return lambda * 1e9; });
  tf(state, weights);

  for (std::size_t z = 0; z < 2; ++z)
  {
    for (std::size_t col = 0; col < 2; ++col)
    {
      EXPECT_DOUBLE_EQ(weights(0, z, col), 225.0);
      EXPECT_DOUBLE_EQ(weights(1, z, col), 275.0);
      EXPECT_DOUBLE_EQ(weights(2, z, col), 350.0);
    }
  }
}

TEST(Factories, WrapAnalyticWavelengthAndTemperature)
{
  auto state = make_state();
  auto weights = make_weights();

  // f(lambda, T) = T  -  verify T is indexed by (z, col)
  auto tf = tuvx::wrap_analytic([](double /*lambda*/, double temp) { return temp; });
  tf(state, weights);

  for (std::size_t wl = 0; wl < 3; ++wl)
  {
    EXPECT_DOUBLE_EQ(weights(wl, 0, 0), 250.0);
    EXPECT_DOUBLE_EQ(weights(wl, 1, 0), 270.0);
    EXPECT_DOUBLE_EQ(weights(wl, 0, 1), 260.0);
    EXPECT_DOUBLE_EQ(weights(wl, 1, 1), 280.0);
  }
}

TEST(Factories, FromData)
{
  auto state = make_state();
  auto weights = make_weights();

  tuvx::Array1D<double> vals(3);
  vals[0] = 1.0e-22;
  vals[1] = 2.0e-22;
  vals[2] = 3.0e-22;

  const auto tf = tuvx::from_data(vals);
  tf(state, weights);

  for (std::size_t z = 0; z < 2; ++z)
  {
    for (std::size_t col = 0; col < 2; ++col)
    {
      EXPECT_DOUBLE_EQ(weights(0, z, col), 1.0e-22);
      EXPECT_DOUBLE_EQ(weights(1, z, col), 2.0e-22);
      EXPECT_DOUBLE_EQ(weights(2, z, col), 3.0e-22);
    }
  }
}

TEST(Factories, TemperatureInterpolation)
{
  auto state = make_state();
  auto weights = make_weights();

  // interp(wl_idx, T) = wl_idx * 10.0 + T / 100.0
  auto tf = tuvx::temperature_interpolation(
      [](std::size_t wl_idx, double temp)
      { return (static_cast<double>(wl_idx) * 10.0) + (temp / 100.0); });
  tf(state, weights);

  for (std::size_t wl = 0; wl < 3; ++wl)
  {
    const double wl_part = static_cast<double>(wl) * 10.0;
    EXPECT_DOUBLE_EQ(weights(wl, 0, 0), wl_part + (250.0 / 100.0));
    EXPECT_DOUBLE_EQ(weights(wl, 1, 0), wl_part + (270.0 / 100.0));
    EXPECT_DOUBLE_EQ(weights(wl, 0, 1), wl_part + (260.0 / 100.0));
    EXPECT_DOUBLE_EQ(weights(wl, 1, 1), wl_part + (280.0 / 100.0));
  }
}

// polynomial_scaling: coeffs[wl x order], computes sum_n coeffs(wl,n) * dT^n
// T_ref = 250 K
//   wl=0: coeffs = [1.0, 0.01, 0.0]  -> P(dT) = 1.0 + 0.01*dT
//   wl=1: coeffs = [2.0, 0.0, 0.001] -> P(dT) = 2.0 + 0.001*dT^2
//   wl=2: coeffs = [0.5, 0.02, 0.0]  -> P(dT) = 0.5 + 0.02*dT
TEST(Factories, PolynomialScaling)
{
  auto state = make_state();
  auto weights = make_weights();

  tuvx::Array2D<double> coeffs(3, 3);
  coeffs(0, 0) = 1.0;
  coeffs(0, 1) = 0.01;
  coeffs(0, 2) = 0.0;
  coeffs(1, 0) = 2.0;
  coeffs(1, 1) = 0.0;
  coeffs(1, 2) = 0.001;
  coeffs(2, 0) = 0.5;
  coeffs(2, 1) = 0.02;
  coeffs(2, 2) = 0.0;

  const auto tf = tuvx::polynomial_scaling(coeffs, 250.0);
  tf(state, weights);

  // col0, z0: T=250 -> dT=0
  EXPECT_DOUBLE_EQ(weights(0, 0, 0), 1.0);
  EXPECT_DOUBLE_EQ(weights(1, 0, 0), 2.0);
  EXPECT_DOUBLE_EQ(weights(2, 0, 0), 0.5);

  // col0, z1: T=270 -> dT=20
  EXPECT_DOUBLE_EQ(weights(0, 1, 0), 1.0 + (0.01 * 20.0));
  EXPECT_DOUBLE_EQ(weights(1, 1, 0), 2.0 + (0.001 * 20.0 * 20.0));
  EXPECT_DOUBLE_EQ(weights(2, 1, 0), 0.5 + (0.02 * 20.0));

  // col1, z0: T=260 -> dT=10
  EXPECT_DOUBLE_EQ(weights(0, 0, 1), 1.0 + (0.01 * 10.0));
  EXPECT_DOUBLE_EQ(weights(1, 0, 1), 2.0 + (0.001 * 10.0 * 10.0));
  EXPECT_DOUBLE_EQ(weights(2, 0, 1), 0.5 + (0.02 * 10.0));
}

TEST(Factories, ExponentialScaling)
{
  auto state = make_state();
  auto weights = make_weights();

  // coeffs: order-1 polynomial in exponent: exp(0.0 + 0.001*dT)
  tuvx::Array2D<double> coeffs(3, 2);
  for (std::size_t wl = 0; wl < 3; ++wl)
  {
    coeffs(wl, 0) = 0.0;
    coeffs(wl, 1) = 0.001;
  }

  const auto tf = tuvx::exponential_scaling(coeffs, 250.0);
  tf(state, weights);

  for (std::size_t wl = 0; wl < 3; ++wl)
  {
    EXPECT_DOUBLE_EQ(weights(wl, 0, 0), std::exp(0.0));
    EXPECT_NEAR(weights(wl, 1, 0), std::exp(0.001 * 20.0), 1e-12);
    EXPECT_NEAR(weights(wl, 0, 1), std::exp(0.001 * 10.0), 1e-12);
  }
}

TEST(Factories, LinearCorrection)
{
  auto state = make_state();
  auto weights = make_weights();

  tuvx::Array1D<double> base(3);
  tuvx::Array1D<double> slope(3);
  base[0] = 1.0e-22;
  base[1] = 2.0e-22;
  base[2] = 3.0e-22;
  slope[0] = 1.0e-24;
  slope[1] = 2.0e-24;
  slope[2] = 3.0e-24;

  const auto tf = tuvx::linear_correction(base, slope, 250.0);
  tf(state, weights);

  // col0, z1: T=270 -> dT=20
  EXPECT_DOUBLE_EQ(weights(0, 1, 0), 1.0e-22 + (1.0e-24 * 20.0));
  EXPECT_DOUBLE_EQ(weights(1, 1, 0), 2.0e-22 + (2.0e-24 * 20.0));
  // col0, z0: T=250 -> dT=0
  EXPECT_DOUBLE_EQ(weights(2, 0, 0), 3.0e-22);
}

// stern_volmer: phi = phi0 / (1 + k*M*phi0)
// phi0=0.4, k=1.0e3 m3/mol, M varies per (z, col)
TEST(Factories, SternVolmer)
{
  auto state = make_state();
  auto weights = make_weights();

  const double phi0 = 0.4;
  const double k_q = 1.0e3;

  const auto tf = tuvx::stern_volmer(phi0, k_q);
  tf(state, weights);

  for (std::size_t wl = 0; wl < 3; ++wl)
  {
    for (std::size_t z = 0; z < 2; ++z)
    {
      for (std::size_t col = 0; col < 2; ++col)
      {
        const double M = state.air_density_(z, col);
        const double expected = phi0 / (1.0 + (k_q * M * phi0));
        EXPECT_NEAR(weights(wl, z, col), expected, 1e-15);
      }
    }
  }
}

TEST(Factories, Parameterized)
{
  auto state = make_state();
  auto weights = make_weights();

  // weight = wl_idx * 100 + z_idx * 10 + col_idx
  auto tf = tuvx::parameterized(
      [](std::size_t wl, std::size_t z, std::size_t col, double /*temp*/, double /*air*/)
      { return static_cast<double>((wl * 100) + (z * 10) + col); });
  tf(state, weights);

  for (std::size_t wl = 0; wl < 3; ++wl)
  {
    for (std::size_t z = 0; z < 2; ++z)
    {
      for (std::size_t col = 0; col < 2; ++col)
      {
        EXPECT_DOUBLE_EQ(weights(wl, z, col), static_cast<double>((wl * 100) + (z * 10) + col));
      }
    }
  }
}
