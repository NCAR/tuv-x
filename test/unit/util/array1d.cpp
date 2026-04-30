// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#include <tuvx/util/array1d.hpp>

#include <gtest/gtest.h>

TEST(Array1D, DefaultConstructor)
{
  tuvx::Array1D<int> a;
  EXPECT_EQ(a.Size(), 0);
}

TEST(Array1D, SizedConstructor)
{
  tuvx::Array1D<int> a(5);
  EXPECT_EQ(a.Size(), 5);
  for (std::size_t i = 0; i < a.Size(); ++i)
  {
    EXPECT_EQ(a[i], 0);
  }
}

TEST(Array1D, FillConstructor)
{
  tuvx::Array1D<double> a(3, 42.0);
  EXPECT_EQ(a.Size(), 3);
  for (std::size_t i = 0; i < a.Size(); ++i)
  {
    EXPECT_DOUBLE_EQ(a[i], 42.0);
  }
}

TEST(Array1D, ElementAccess)
{
  tuvx::Array1D<int> a(4);
  a[0] = 10;
  a[1] = 20;
  a[2] = 30;
  a[3] = 40;
  EXPECT_EQ(a[0], 10);
  EXPECT_EQ(a[1], 20);
  EXPECT_EQ(a[2], 30);
  EXPECT_EQ(a[3], 40);
}

TEST(Array1D, ConstAccess)
{
  tuvx::Array1D<int> a(3);
  a[0] = 1;
  a[1] = 2;
  a[2] = 3;
  const auto &ca = a;
  EXPECT_EQ(ca[0], 1);
  EXPECT_EQ(ca[1], 2);
  EXPECT_EQ(ca[2], 3);
}

TEST(Array1D, Iteration)
{
  tuvx::Array1D<int> a(5);
  int count = 0;
  for (auto &x : a)
  {
    x = count++;
  }
  EXPECT_EQ(a[0], 0);
  EXPECT_EQ(a[1], 1);
  EXPECT_EQ(a[2], 2);
  EXPECT_EQ(a[3], 3);
  EXPECT_EQ(a[4], 4);
}

TEST(Array1D, ConstIteration)
{
  tuvx::Array1D<int> a(3, 7);
  const auto &ca = a;
  int sum = 0;
  for (const auto &x : ca)
  {
    sum += x;
  }
  EXPECT_EQ(sum, 21);
}

TEST(Array1D, AsVector)
{
  tuvx::Array1D<int> a(3);
  a[0] = 10;
  a[1] = 20;
  a[2] = 30;
  auto &v = a.AsVector();
  EXPECT_EQ(v.size(), 3);
  EXPECT_EQ(v[0], 10);
  EXPECT_EQ(v[1], 20);
  EXPECT_EQ(v[2], 30);

  const auto &ca = a;
  const auto &cv = ca.AsVector();
  EXPECT_EQ(cv.size(), 3);
}
