// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
//
// Regression tests for the analytic spectral weights in the example fixed
// configuration. Each test loads the reference CSV (see reference/README for
// provenance) and compares against the C++ output.
//
// Spectral weights are unitless and wavelength-only, so the reference CSVs have
// two columns: wavelength_m, weight. Tolerance: relative 1e-10.
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
    double weight_{};
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
      ss >> row.wavelength_m_ >> comma >> row.weight_;
      rows.push_back(row);
    }
    return rows;
  }

  std::string RefDir()
  {
    const char *env = std::getenv("TUVX_REGRESSION_REFERENCE_DIR");
    return std::string(env != nullptr ? env : TUVX_REGRESSION_REFERENCE_DIR) + "/";
  }

  void CheckAgainstRef(const std::string &csv_path, const tuvx::TransformFunc<> &tf)
  {
    const auto ref = LoadCsv(csv_path);
    const auto n = ref.size();

    tuvx::AtmosphericState<> state;
    state.wavelength_grid_ = tuvx::Grid<tuvx::Array2D<double>>("m", 1, n);
    for (std::size_t i = 0; i < n; ++i)
    {
      state.wavelength_grid_.mid_points_(i, 0) = ref[i].wavelength_m_;
    }
    state.height_grid_ = tuvx::Grid<tuvx::Array2D<double>>("m", 1, 1);
    state.temperature_ = tuvx::Array2D<double>(1, 1);
    state.pressure_ = tuvx::Array2D<double>(1, 1);
    state.air_density_ = tuvx::Array2D<double>(1, 1);

    tuvx::Array3D<double> weights(n, 1, 1);
    tf(state, weights);

    for (std::size_t i = 0; i < n; ++i)
    {
      const double got = weights(i, 0, 0);
      const double expected = ref[i].weight_;
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

TEST(SpectralWeights, StandardHumanErythema)
{
  CheckAgainstRef(
      RefDir() + "standard_human_erythema.csv", tuvx::fixed_configuration::spectral_weights::standard_human_erythema());
}

TEST(SpectralWeights, UVIndex)
{
  CheckAgainstRef(RefDir() + "uv_index.csv", tuvx::fixed_configuration::spectral_weights::uv_index());
}

TEST(SpectralWeights, ScupMice)
{
  CheckAgainstRef(RefDir() + "scup_mice.csv", tuvx::fixed_configuration::spectral_weights::scup_mice());
}

TEST(SpectralWeights, ExpDecay)
{
  CheckAgainstRef(RefDir() + "exp_decay.csv", tuvx::fixed_configuration::spectral_weights::exp_decay());
}

TEST(SpectralWeights, PAR)
{
  CheckAgainstRef(RefDir() + "par.csv", tuvx::fixed_configuration::spectral_weights::par());
}

TEST(SpectralWeights, Gaussian)
{
  CheckAgainstRef(RefDir() + "gaussian.csv", tuvx::fixed_configuration::spectral_weights::gaussian(305.0));
}

TEST(SpectralWeights, PlantDamage)
{
  CheckAgainstRef(RefDir() + "plant_damage.csv", tuvx::fixed_configuration::spectral_weights::plant_damage());
}

TEST(SpectralWeights, PlantDamageFlintCaldwell)
{
  CheckAgainstRef(
      RefDir() + "plant_damage_flint_caldwell.csv",
      tuvx::fixed_configuration::spectral_weights::plant_damage_flint_caldwell());
}

TEST(SpectralWeights, PlantDamageFlintCaldwellExt)
{
  CheckAgainstRef(
      RefDir() + "plant_damage_flint_caldwell_ext.csv",
      tuvx::fixed_configuration::spectral_weights::plant_damage_flint_caldwell_ext());
}

TEST(SpectralWeights, PhytoplanktonBoucher)
{
  CheckAgainstRef(
      RefDir() + "phytoplankton_boucher.csv", tuvx::fixed_configuration::spectral_weights::phytoplankton_boucher());
}
