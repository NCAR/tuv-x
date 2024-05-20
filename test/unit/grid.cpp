// Copyright (C) 2023-2024 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
//
// Tests for tuvx::Grid.
#include <gtest/gtest.h>

#include <tuvx/grid.hpp>

TEST(Grid, Constructor) {
  tuvx::Grid<tuvx::Array2D<int>> a("foos", 2, 3);
  EXPECT_EQ(a.units(), "foos");
  EXPECT_EQ(a.number_of_columns(), 2);
  EXPECT_EQ(a.number_of_cells(), 3);
  EXPECT_EQ(a.number_of_edges(), 4);
  EXPECT_FALSE(a.is_constant());
  EXPECT_EQ(a.mid_points_.Size1(), 3);
  EXPECT_EQ(a.mid_points_.Size2(), 2);
  EXPECT_EQ(a.edges_.Size1(), 4);
  EXPECT_EQ(a.edges_.Size2(), 2);
  
  tuvx::Grid<tuvx::Array2D<int>> b("bars", 4);
  EXPECT_STREQ(b.units().c_str(), "bars");
  EXPECT_EQ(b.number_of_columns(), 1);
  EXPECT_EQ(b.number_of_cells(), 4);
  EXPECT_EQ(b.number_of_edges(), 5);
  EXPECT_TRUE(b.is_constant());
  EXPECT_EQ(b.mid_points_.Size1(), 4);
  EXPECT_EQ(b.mid_points_.Size2(), 1);
  EXPECT_EQ(b.edges_.Size1(), 5);
  EXPECT_EQ(b.edges_.Size2(), 1);
}