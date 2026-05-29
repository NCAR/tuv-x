// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#include <tuvx/radiative_transfer/spherical_geometry.hpp>

#include <gtest/gtest.h>

#include <cmath>
#include <optional>
#include <vector>

// For overhead sun (SZA = 0), dsdh(level, layer) is analytically exactly 1.0
// for all layers reached from a given level.  This follows because the slant
// path through a spherical shell equals the vertical depth when the beam is
// purely vertical.
TEST(SphericalGeometry, OverheadSunDsdh)
{
  // 3-layer atmosphere: ground at 0 m, mid at 5 km, 10 km, TOA at 15 km
  std::vector<double> edges = { 0.0, 5000.0, 10000.0, 15000.0 };
  tuvx::SphericalGeometry geom;
  geom.SetParameters(0.0, edges);  // SZA = 0 rad

  const std::size_t N_LAYERS = edges.size() - 1;
  // At SZA=0, nid[i] = i for all levels
  for (std::size_t i = 0; i <= N_LAYERS; ++i)
  {
    EXPECT_EQ(geom.nid_[i], std::optional<std::size_t>(i)) << "level " << i;
  }
  // dsdh[i][j] = 1.0 for j < i (layers above level i)
  for (std::size_t i = 1; i <= N_LAYERS; ++i)
  {
    for (std::size_t j = 0; j < i; ++j)
    {
      EXPECT_NEAR(geom.dsdh_[i][j], 1.0, 1.0e-10) << "level " << i << " layer " << j;
    }
  }
}

// At SZA = 0 and constant optical depth per layer, the slant optical depth
// at level i equals the sum of taun[0..i-1] (purely vertical path).
TEST(SphericalGeometry, SlantOpticalDepthOverheadSun)
{
  std::vector<double> edges = { 0.0, 1000.0, 2000.0, 3000.0 };
  tuvx::SphericalGeometry geom;
  geom.SetParameters(0.0, edges);

  std::vector<double> taun = { 0.1, 0.2, 0.3 };  // top-to-bottom layers

  // At level 0 (TOA), path = 0
  EXPECT_NEAR(tuvx::SlantOpticalDepth(0, geom.nid_[0], geom.dsdh_[0], taun), 0.0, 1.0e-12);
  // Level 1: sum of layer 0 = 0.1
  EXPECT_NEAR(tuvx::SlantOpticalDepth(1, geom.nid_[1], geom.dsdh_[1], taun), 0.1, 1.0e-10);
  // Level 2: layers 0+1 = 0.3
  EXPECT_NEAR(tuvx::SlantOpticalDepth(2, geom.nid_[2], geom.dsdh_[2], taun), 0.3, 1.0e-10);
  // Level 3 (ground): layers 0+1+2 = 0.6
  EXPECT_NEAR(tuvx::SlantOpticalDepth(3, geom.nid_[3], geom.dsdh_[3], taun), 0.6, 1.0e-10);
}

// nullopt nid (below tangent height) returns infinity.
TEST(SphericalGeometry, SlantOpticalDepthBelowTangentHeight)
{
  std::vector<double> taun = { 0.1 };
  std::vector<double> slpath = { 1.0 };
  double result = tuvx::SlantOpticalDepth(0, std::nullopt, slpath, taun);
  EXPECT_TRUE(std::isinf(result));
}

// For overhead sun, slant_optical_depth at level 0 with nid=0 should be 0.
TEST(SphericalGeometry, SlantOpticalDepthAtToa)
{
  std::vector<double> taun = { 0.5, 0.3 };
  std::vector<double> slpath = { 1.0, 1.0 };
  // nid[0] = 0 for overhead sun → no layers crossed → path = 0
  EXPECT_DOUBLE_EQ(tuvx::SlantOpticalDepth(0, std::size_t{ 0 }, slpath, taun), 0.0);
}
