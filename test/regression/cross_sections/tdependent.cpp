// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
//
// Regression tests for the temperature-dependent cross-sections in the example
// fixed configuration (N2O, Cl2, H2O2, CHBr3). Each reference CSV carries a
// temperature axis: rows are (wavelength_m, temperature_K, cross_section_m2).
//
// The state is built with the distinct wavelengths on the wavelength axis and
// the distinct temperatures on the column axis (one height layer), so a single
// transform evaluation covers every (wavelength, temperature) pair in the file.
//
// Tolerance: relative 1e-10 -- both sides evaluate the same IEEE 754 formulas.
#include <tuvx/fixed_configuration.hpp>

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{
  struct RefRow
  {
    double wavelength_m_{};
    double temperature_k_{};
    double cross_section_m2_{};
  };

  std::vector<RefRow> LoadCsv(const std::string &path)
  {
    std::ifstream f(path);
    if (!f)
    {
      throw std::runtime_error("Cannot open reference file: " + path);
    }
    std::string line;
    std::getline(f, line);  // skip header
    std::vector<RefRow> rows;
    while (std::getline(f, line))
    {
      if (line.empty())
      {
        continue;
      }
      std::istringstream ss(line);
      RefRow row{};
      char comma = 0;
      ss >> row.wavelength_m_ >> comma >> row.temperature_k_ >> comma >> row.cross_section_m2_;
      rows.push_back(row);
    }
    return rows;
  }

  std::string RefDir()
  {
    const char *env = std::getenv("TUVX_REGRESSION_REFERENCE_DIR");
    return std::string(env != nullptr ? env : TUVX_REGRESSION_REFERENCE_DIR) + "/";
  }

  // Distinct values, preserving first-seen order.
  std::vector<double> Distinct(const std::vector<RefRow> &rows, double RefRow::*field)
  {
    std::vector<double> out;
    for (const auto &row : rows)
    {
      const double v = row.*field;
      if (std::ranges::find(out, v) == out.end())
      {
        out.push_back(v);
      }
    }
    return out;
  }

  std::size_t IndexOf(const std::vector<double> &values, double v)
  {
    return static_cast<std::size_t>(std::ranges::find(values, v) - values.begin());
  }

  void CheckAgainstRef(const std::string &csv_path, const tuvx::TransformFunc<> &tf, double rel_tol = 1.0e-10)
  {
    const auto ref = LoadCsv(csv_path);
    const auto wavelengths = Distinct(ref, &RefRow::wavelength_m_);
    const auto temperatures = Distinct(ref, &RefRow::temperature_k_);
    const auto n_wl = wavelengths.size();
    const auto n_col = temperatures.size();

    // wavelength axis -> rows; temperature axis -> columns; single height layer.
    tuvx::AtmosphericState<> state;
    state.wavelength_grid_ = tuvx::Grid<tuvx::Array2D<double>>("m", 1, n_wl);
    for (std::size_t i = 0; i < n_wl; ++i)
    {
      state.wavelength_grid_.mid_points_(i, 0) = wavelengths[i];
    }
    state.height_grid_ = tuvx::Grid<tuvx::Array2D<double>>("m", 1, n_col);
    state.temperature_ = tuvx::Array2D<double>(1, n_col);
    state.pressure_ = tuvx::Array2D<double>(1, n_col);
    state.air_density_ = tuvx::Array2D<double>(1, n_col);
    for (std::size_t col = 0; col < n_col; ++col)
    {
      state.temperature_(0, col) = temperatures[col];
    }

    tuvx::Array3D<double> weights(n_wl, 1, n_col);
    tf(state, weights);

    for (const auto &row : ref)
    {
      const std::size_t wl = IndexOf(wavelengths, row.wavelength_m_);
      const std::size_t col = IndexOf(temperatures, row.temperature_k_);
      const double got = weights(wl, 0, col);
      const double expected = row.cross_section_m2_;
      if (expected == 0.0)
      {
        EXPECT_DOUBLE_EQ(got, 0.0) << "at wavelength " << row.wavelength_m_ << " m, T " << row.temperature_k_ << " K";
      }
      else
      {
        const double rel_err = std::abs((got - expected) / expected);
        EXPECT_LT(rel_err, rel_tol) << "at wavelength " << row.wavelength_m_ << " m, T " << row.temperature_k_
                                    << " K  got=" << got << "  expected=" << expected;
      }
    }
  }
}  // namespace

TEST(TemperatureDependentCrossSections, N2O)
{
  CheckAgainstRef(RefDir() + "n2o.csv", tuvx::fixed_configuration::n2o());
}

TEST(TemperatureDependentCrossSections, Cl2)
{
  CheckAgainstRef(RefDir() + "cl2.csv", tuvx::fixed_configuration::cl2());
}

TEST(TemperatureDependentCrossSections, H2O2)
{
  // Looser tolerance for H2O2 only: its degree-7 polynomial in wavelength has
  // large alternating coefficients and the terms nearly cancel, so the result
  // is ill-conditioned. Different platforms' libm / FP contraction evaluate it
  // slightly differently (worst observed ~1.4e-10 relative, deterministic per
  // platform). 1e-8 still asserts ~8 significant figures of agreement.
  CheckAgainstRef(RefDir() + "h2o2.csv", tuvx::fixed_configuration::h2o2(), 1.0e-8);
}

TEST(TemperatureDependentCrossSections, CHBr3)
{
  CheckAgainstRef(RefDir() + "chbr3.csv", tuvx::fixed_configuration::chbr3());
}
