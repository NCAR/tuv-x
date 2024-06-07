/* Copyright (C) 2023-2024 National Center for Atmospheric Research
 *
 * SPDX-License-Identifier: Apache-2.0
 */
#include <tuvx/grid.hpp>

#include <gtest/gtest.h>

TEST(Grid, Constructor)
{
  tuvx::Grid<tuvx::Array2D<int>> a("foos", 2, 3);
  EXPECT_EQ(a.Units(), "foos");
  EXPECT_EQ(a.NumberOfColumns(), 2);
  EXPECT_EQ(a.NumberOfSections(), 3);
  EXPECT_EQ(a.NumberOfEdges(), 4);
  EXPECT_FALSE(a.IsConstant());
  EXPECT_EQ(a.mid_points_.Size1(), 3);
  EXPECT_EQ(a.mid_points_.Size2(), 2);
  EXPECT_EQ(a.edges_.Size1(), 4);
  EXPECT_EQ(a.edges_.Size2(), 2);

  tuvx::Grid<tuvx::Array2D<int>> b("bars", 4);
  EXPECT_STREQ(b.Units().c_str(), "bars");
  EXPECT_EQ(b.NumberOfColumns(), 1);
  EXPECT_EQ(b.NumberOfSections(), 4);
  EXPECT_EQ(b.NumberOfEdges(), 5);
  EXPECT_TRUE(b.IsConstant());
  EXPECT_EQ(b.mid_points_.Size1(), 4);
  EXPECT_EQ(b.mid_points_.Size2(), 1);
  EXPECT_EQ(b.edges_.Size1(), 5);
  EXPECT_EQ(b.edges_.Size2(), 1);
}