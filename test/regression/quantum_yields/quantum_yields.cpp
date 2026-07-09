// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
//
// Regression tests for the pure-wavelength quantum yields in the example fixed
// configuration. Each test loads the reference CSV (see reference/README) and
// compares against the C++ output.
//
// Quantum yields are unitless and wavelength-only here, so the reference CSVs
// have two columns: wavelength_m, quantum_yield. Tolerance: relative 1e-10.
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
    double quantum_yield_{};
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
      ss >> row.wavelength_m_ >> comma >> row.quantum_yield_;
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
      const double expected = ref[i].quantum_yield_;
      if (expected == 0.0)
      {
        EXPECT_DOUBLE_EQ(got, 0.0) << "at wavelength " << ref[i].wavelength_m_ << " m";
      }
      else
      {
        const double rel_err = std::abs((got - expected) / expected);
        EXPECT_LT(rel_err, 1.0e-10)
            << "at wavelength " << ref[i].wavelength_m_ << " m  got=" << got << "  expected=" << expected;
      }
    }
  }
}  // namespace

TEST(QuantumYields, CloClO1D)
{
  CheckAgainstRef(RefDir() + "clo_cl_o1d.csv", tuvx::fixed_configuration::quantum_yields::clo_cl_o1d());
}

TEST(QuantumYields, CloClO3P)
{
  CheckAgainstRef(RefDir() + "clo_cl_o3p.csv", tuvx::fixed_configuration::quantum_yields::clo_cl_o3p());
}

TEST(QuantumYields, ClONO2ClNO3)
{
  CheckAgainstRef(RefDir() + "clono2_cl_no3.csv", tuvx::fixed_configuration::quantum_yields::clono2_cl_no3());
}

TEST(QuantumYields, ClONO2ClONO2)
{
  CheckAgainstRef(RefDir() + "clono2_clo_no2.csv", tuvx::fixed_configuration::quantum_yields::clono2_clo_no2());
}

TEST(QuantumYields, HO2OHO)
{
  CheckAgainstRef(RefDir() + "ho2_oh_o.csv", tuvx::fixed_configuration::quantum_yields::ho2_oh_o());
}
