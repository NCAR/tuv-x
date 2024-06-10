// Copyright (C) 2023-2024 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#include <tuvx/util/array3d.hpp>

#include <gtest/gtest.h>

TEST(Array3D, Constructor)
{
  tuvx::Array3D<int> a(2, 3, 4);
  EXPECT_EQ(a.Size1(), 2);
  EXPECT_EQ(a.Size2(), 3);
  EXPECT_EQ(a.Size3(), 4);
  int count = 10;
  for (auto &x : a)
  {
    x = count++;
  }
  a(1, 0, 0) = 45;
  a(1, 2, 3) = 100;
  EXPECT_EQ(a(0, 0, 0), 10);
  EXPECT_EQ(a(0, 0, 1), 11);
  EXPECT_EQ(a(0, 0, 2), 12);
  EXPECT_EQ(a(0, 0, 3), 13);
  EXPECT_EQ(a(0, 1, 0), 14);
  EXPECT_EQ(a(0, 1, 1), 15);
  EXPECT_EQ(a(0, 1, 2), 16);
  EXPECT_EQ(a(0, 1, 3), 17);
  EXPECT_EQ(a(0, 2, 0), 18);
  EXPECT_EQ(a(0, 2, 1), 19);
  EXPECT_EQ(a(0, 2, 2), 20);
  EXPECT_EQ(a(0, 2, 3), 21);
  EXPECT_EQ(a(1, 0, 0), 45);
  EXPECT_EQ(a(1, 0, 1), 23);
  EXPECT_EQ(a(1, 0, 2), 24);
  EXPECT_EQ(a(1, 0, 3), 25);
  EXPECT_EQ(a(1, 1, 0), 26);
  EXPECT_EQ(a(1, 1, 1), 27);
  EXPECT_EQ(a(1, 1, 2), 28);
  EXPECT_EQ(a(1, 1, 3), 29);
  EXPECT_EQ(a(1, 2, 0), 30);
  EXPECT_EQ(a(1, 2, 1), 31);
  EXPECT_EQ(a(1, 2, 2), 32);
  EXPECT_EQ(a(1, 2, 3), 100);
}