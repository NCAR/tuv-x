// Copyright (C) 2023-2024 National Science Foundation-National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#include <tuvx/util/array2d.hpp>

#include <gtest/gtest.h>

TEST(Array2D, Constructor)
{
  tuvx::Array2D<int> a(2, 3);
  EXPECT_EQ(a.Size1(), 2);
  EXPECT_EQ(a.Size2(), 3);
  int count = 10;
  for (auto &x : a)
  {
    x = count++;
  }
  a(1, 0) = 45;
  a(1, 2) = 100;
  EXPECT_EQ(a(0, 0), 10);
  EXPECT_EQ(a(0, 1), 11);
  EXPECT_EQ(a(0, 2), 12);
  EXPECT_EQ(a(1, 0), 45);
  EXPECT_EQ(a(1, 1), 14);
  EXPECT_EQ(a(1, 2), 100);
}