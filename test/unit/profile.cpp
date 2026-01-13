// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#include <tuvx/profile.hpp>

#include <gtest/gtest.h>

TEST(Profile, Constructor)
{
  tuvx::Grid<tuvx::Array2D<int>> grid("foos", 2, 3);
  tuvx::Profile<tuvx::Array2D<int>> a("bars", 2, grid);
  EXPECT_EQ(a.Units(), "bars");
  EXPECT_EQ(a.mid_point_values_.Size1(), 3);
  EXPECT_EQ(a.mid_point_values_.Size2(), 2);
  EXPECT_EQ(a.edge_values_.Size1(), 4);
  EXPECT_EQ(a.edge_values_.Size2(), 2);
}