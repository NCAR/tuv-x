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
#include <limits>
#include <optional>
#include <sstream>
#include <string>
#include <vector>

// ── CSV reader ────────────────────────────────────────────────────────────────

struct ReferenceLevel
{
  std::size_t              level;
  std::optional<std::size_t> nid;          // nullopt if CSV contains -1
  std::vector<double>      dsdh;
  double                   slant_od;       // Inf if CSV contains "Inf"
};

static std::vector<ReferenceLevel> LoadReference(const std::string& path, std::size_t n_layers)
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
    row.dsdh.resize(n_layers);

    std::getline(ss, token, ',');
    row.level = std::stoul(token);

    std::getline(ss, token, ',');
    const int nid_int = std::stoi(token);
    row.nid = (nid_int < 0) ? std::nullopt : std::optional<std::size_t>(static_cast<std::size_t>(nid_int));

    for (std::size_t j = 0; j < n_layers; ++j)
    {
      std::getline(ss, token, ',');
      row.dsdh[j] = std::stod(token);
    }

    std::getline(ss, token, ',');
    row.slant_od = (token.find("Inf") != std::string::npos)
                       ? std::numeric_limits<double>::infinity()
                       : std::stod(token);

    rows.push_back(row);
  }
  return rows;
}

// ── Helper ────────────────────────────────────────────────────────────────────

static std::string ReferenceFile(double sza_rad)
{
  const char* src_dir = std::getenv("TUVX_REGRESSION_REFERENCE_DIR");
  std::string dir     = src_dir ? src_dir : TUVX_REGRESSION_REFERENCE_DIR;
  char        buf[64];
  std::snprintf(buf, sizeof(buf), "spherical_geometry_sza_%.4f.csv", sza_rad);
  return dir + "/" + buf;
}

// ── Test runner ───────────────────────────────────────────────────────────────

static void RunRegressionTest(double sza_rad, double dsdh_tol, double slant_od_tol)
{
  static constexpr std::size_t n_layers = 3;

  const std::vector<double> alt_edges = { 0.0, 1000.0, 2000.0, 3000.0 };
  const std::vector<double> taun      = { 0.1, 0.2, 0.3 };

  tuvx::SphericalGeometry geom;
  geom.SetParameters(sza_rad, alt_edges);

  const std::string path      = ReferenceFile(sza_rad);
  const auto        reference = LoadReference(path, n_layers);
  ASSERT_EQ(reference.size(), n_layers + 1) << "Expected " << n_layers + 1 << " rows in " << path;

  for (const auto& ref : reference)
  {
    const std::size_t lev = ref.level;
    EXPECT_EQ(geom.nid_[lev], ref.nid) << "nid mismatch at level " << lev << " sza=" << sza_rad;

    for (std::size_t j = 0; j < n_layers; ++j)
    {
      EXPECT_NEAR(geom.dsdh_[lev][j], ref.dsdh[j], dsdh_tol)
          << "dsdh[" << lev << "][" << j << "] mismatch sza=" << sza_rad;
    }

    const double slant_od = tuvx::SlantOpticalDepth(lev, geom.nid_[lev], geom.dsdh_[lev], taun);
    if (std::isinf(ref.slant_od))
    {
      EXPECT_TRUE(std::isinf(slant_od)) << "expected Inf at level " << lev << " sza=" << sza_rad;
    }
    else
    {
      EXPECT_NEAR(slant_od, ref.slant_od, slant_od_tol)
          << "slant_od mismatch at level " << lev << " sza=" << sza_rad;
    }
  }
}

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
