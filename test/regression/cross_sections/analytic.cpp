// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
//
// Regression tests for the analytic cross-sections in the example fixed
// configuration.  Each test loads the reference CSV (see reference/README for
// provenance) and compares against the C++ output assembled from the general
// analytic transform forms.
//
// Tolerance: relative 1e-10 -- both sides evaluate the same IEEE 754 formulas.
#include <tuvx/fixed_configuration.hpp>

#include <gtest/gtest.h>

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
      ss >> row.wavelength_m_ >> comma >> row.cross_section_m2_;
      rows.push_back(row);
    }
    return rows;
  }

  std::string RefDir()
  {
    const char *env = std::getenv("TUVX_REGRESSION_REFERENCE_DIR");
    return std::string(env != nullptr ? env : TUVX_REGRESSION_REFERENCE_DIR) + "/";
  }

  tuvx::AtmosphericState<> MakeState(const std::vector<RefRow> &ref)
  {
    tuvx::AtmosphericState<> state;
    const auto n = ref.size();
    state.wavelength_grid_ = tuvx::Grid<tuvx::Array2D<double>>("m", 1, n);
    for (std::size_t i = 0; i < n; ++i)
    {
      state.wavelength_grid_.mid_points_(i, 0) = ref[i].wavelength_m_;
    }
    state.wavelength_grid_.edges_(0, 0) = ref[0].wavelength_m_ - 5.0e-9;
    for (std::size_t i = 0; i < n; ++i)
    {
      state.wavelength_grid_.edges_(i + 1, 0) = ref[i].wavelength_m_ + 5.0e-9;
    }
    state.height_grid_ = tuvx::Grid<tuvx::Array2D<double>>("m", 1, 1);
    state.temperature_ = tuvx::Array2D<double>(1, 1);
    state.pressure_ = tuvx::Array2D<double>(1, 1);
    state.air_density_ = tuvx::Array2D<double>(1, 1);
    return state;
  }

  void CheckAgainstRef(const std::string &csv_path, const tuvx::TransformFunc<> &tf)
  {
    const auto ref = LoadCsv(csv_path);
    auto state = MakeState(ref);
    const auto n = ref.size();
    tuvx::Array3D<double> weights(n, 1, 1);

    tf(state, weights);

    for (std::size_t i = 0; i < n; ++i)
    {
      const double got = weights(i, 0, 0);
      const double expected = ref[i].cross_section_m2_;
      if (expected == 0.0)
      {
        EXPECT_DOUBLE_EQ(got, 0.0) << "at wavelength " << ref[i].wavelength_m_ << " m";
      }
      else
      {
        const double rel_err = std::abs((got - expected) / expected);
        EXPECT_LT(rel_err, 1.0e-10) << "at wavelength " << ref[i].wavelength_m_ << " m  got=" << got
                                    << "  expected=" << expected;
      }
    }
  }
}  // namespace

TEST(AnalyticCrossSections, Rayleigh)
{
  CheckAgainstRef(RefDir() + "rayleigh.csv", tuvx::fixed_configuration::rayleigh());
}

TEST(AnalyticCrossSections, HOBr)
{
  CheckAgainstRef(RefDir() + "hobr.csv", tuvx::fixed_configuration::hobr());
}

TEST(AnalyticCrossSections, TButylNitrate)
{
  CheckAgainstRef(RefDir() + "t_butyl_nitrate.csv", tuvx::fixed_configuration::t_butyl_nitrate());
}

TEST(AnalyticCrossSections, NitroxyAcetone)
{
  CheckAgainstRef(RefDir() + "nitroxy_acetone.csv", tuvx::fixed_configuration::nitroxy_acetone());
}

TEST(AnalyticCrossSections, NitroxyEthanol)
{
  CheckAgainstRef(RefDir() + "nitroxy_ethanol.csv", tuvx::fixed_configuration::nitroxy_ethanol());
}
