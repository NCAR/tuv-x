// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
//
// Regression tests for SphericalGeometry::SetParameters() and
// SlantOpticalDepth() against reference CSV data.
//
// Reference CSVs were generated from the C++ implementation at the commit
// recorded in README.md. See README.md for the test atmosphere and
// regeneration steps.
//
// Test atmosphere: 3-layer, altitude edges 0/1/2/3 km, optical depths
//                 0.1/0.2/0.3 (top-to-bottom) for slant_od columns.
#include <tuvx/radiative_transfer/spherical_geometry.hpp>

#include <gtest/gtest.h>

#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <limits>
#include <optional>
#include <sstream>
#include <string>
#include <vector>

// ── CSV reader ────────────────────────────────────────────────────────────────

namespace
{

  struct ReferenceLevel
  {
    std::size_t level_ = 0;
    std::optional<std::size_t> nid_ = std::nullopt;
    std::vector<double> dsdh_ = {};
    double slant_od_ = 0.0;
  };

  std::vector<ReferenceLevel> LoadReference(const std::string& path, std::size_t n_layers)
  {
    std::ifstream file(path);
    EXPECT_TRUE(file.is_open()) << "Cannot open reference file: " << path;

    std::vector<ReferenceLevel> rows;
    std::string line;
    bool header_skipped = false;
    while (std::getline(file, line))
    {
      if (line.empty() || line[0] == '#')
      {
        continue;
      }
      if (!header_skipped)
      {
        header_skipped = true;
        continue;
      }
      std::istringstream ss(line);
      std::string token;
      ReferenceLevel row{};
      row.dsdh_.resize(n_layers);

      std::getline(ss, token, ',');
      row.level_ = std::stoul(token);

      std::getline(ss, token, ',');
      const int NID_INT = std::stoi(token);
      row.nid_ = (NID_INT < 0) ? std::nullopt : std::optional<std::size_t>(static_cast<std::size_t>(NID_INT));

      for (std::size_t j = 0; j < n_layers; ++j)
      {
        std::getline(ss, token, ',');
        row.dsdh_[j] = std::stod(token);
      }

      std::getline(ss, token, ',');
      row.slant_od_ = (token.find("Inf") != std::string::npos) ? std::numeric_limits<double>::infinity() : std::stod(token);

      rows.push_back(row);
    }
    return rows;
  }

  // ── Helpers ───────────────────────────────────────────────────────────────────

  std::string ReferenceFile(double sza_rad)
  {
    const char* src_dir = std::getenv("TUVX_REGRESSION_REFERENCE_DIR");
    const std::string DIR = (src_dir != nullptr) ? src_dir : TUVX_REGRESSION_REFERENCE_DIR;
    std::ostringstream oss;
    oss << std::fixed << std::setprecision(4) << sza_rad;
    return DIR + "/spherical_geometry_sza_" + oss.str() + ".csv";
  }

  void CheckReferenceLevel(
      const tuvx::SphericalGeometry& geom,
      const ReferenceLevel& ref,
      const std::vector<double>& taun,
      std::size_t n_layers,
      double dsdh_tol,
      double slant_od_tol,
      double sza_rad)
  {
    const std::size_t LEV = ref.level_;
    EXPECT_EQ(geom.nid_[LEV], ref.nid_) << "nid mismatch at level " << LEV << " sza=" << sza_rad;

    for (std::size_t j = 0; j < n_layers; ++j)
    {
      EXPECT_NEAR(geom.dsdh_[LEV][j], ref.dsdh_[j], dsdh_tol) << "dsdh[" << LEV << "][" << j << "] mismatch sza=" << sza_rad;
    }

    const double SLANT_OD = tuvx::SlantOpticalDepth(LEV, geom.nid_[LEV], geom.dsdh_[LEV], taun);
    if (std::isinf(ref.slant_od_))
    {
      EXPECT_TRUE(std::isinf(SLANT_OD)) << "expected Inf at level " << LEV << " sza=" << sza_rad;
    }
    else
    {
      EXPECT_NEAR(SLANT_OD, ref.slant_od_, slant_od_tol) << "slant_od mismatch at level " << LEV << " sza=" << sza_rad;
    }
  }

  // ── Test runner ───────────────────────────────────────────────────────────────

  void RunRegressionTest(double sza_rad, double dsdh_tol, double slant_od_tol)
  {
    static constexpr std::size_t N_LAYERS = 3;

    const std::vector<double> ALT_EDGES = { 0.0, 1000.0, 2000.0, 3000.0 };
    const std::vector<double> TAUN = { 0.1, 0.2, 0.3 };

    tuvx::SphericalGeometry geom;
    geom.SetParameters(sza_rad, ALT_EDGES);

    const std::string PATH = ReferenceFile(sza_rad);
    const auto REFERENCE = LoadReference(PATH, N_LAYERS);
    ASSERT_EQ(REFERENCE.size(), N_LAYERS + 1) << "Expected " << N_LAYERS + 1 << " rows in " << PATH;

    for (const auto& ref : REFERENCE)
    {
      CheckReferenceLevel(geom, ref, TAUN, N_LAYERS, dsdh_tol, slant_od_tol, sza_rad);
    }
  }

}  // namespace

// ── Tests ─────────────────────────────────────────────────────────────────────

TEST(SphericalGeometryRegression, Sza0_0)
{
  RunRegressionTest(0.0, 1.0e-12, 1.0e-12);
}

TEST(SphericalGeometryRegression, Sza0_5)
{
  RunRegressionTest(0.5, 1.0e-10, 1.0e-10);
}

TEST(SphericalGeometryRegression, Sza1_0)
{
  RunRegressionTest(1.0, 1.0e-10, 1.0e-10);
}

TEST(SphericalGeometryRegression, Sza1_5)
{
  RunRegressionTest(1.5, 1.0e-10, 1.0e-10);
}

TEST(SphericalGeometryRegression, Sza1_6)
{
  RunRegressionTest(1.6, 1.0e-10, 1.0e-10);
}
