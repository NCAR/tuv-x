// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#include <tuvx/util/array3d.hpp>

#include <tuvx/util/array2d.hpp>

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

TEST(Array3D, ValueType)
{
  static_assert(std::is_same_v<tuvx::Array3D<double>::value_type, double>);
  static_assert(std::is_same_v<tuvx::Array3D<float>::value_type, float>);
}

// Layout [2 × 3 × 2] (wl=2, height=3, col=2)
//   col 0: a(i,j,0), col 1: a(i,j,1)
TEST(Array3D, ColumnView)
{
  tuvx::Array3D<double> a(2, 3, 2);
  // Fill col 0 with sequential values 1..6, col 1 with 10..15
  int c0 = 1;
  int c1 = 10;
  for (std::size_t i = 0; i < 2; ++i)
  {
    for (std::size_t j = 0; j < 3; ++j)
    {
      a(i, j, 0) = static_cast<double>(c0++);
      a(i, j, 1) = static_cast<double>(c1++);
    }
  }

  auto col0 = a.GetColumnView(0);
  EXPECT_DOUBLE_EQ(col0(0, 0), 1.0);
  EXPECT_DOUBLE_EQ(col0(0, 1), 2.0);
  EXPECT_DOUBLE_EQ(col0(0, 2), 3.0);
  EXPECT_DOUBLE_EQ(col0(1, 0), 4.0);
  EXPECT_DOUBLE_EQ(col0(1, 1), 5.0);
  EXPECT_DOUBLE_EQ(col0(1, 2), 6.0);

  // Mutate through view
  col0(1, 0) = 99.0;
  EXPECT_DOUBLE_EQ(a(1, 0, 0), 99.0);

  const auto &ca = a;
  auto cc = ca.GetConstColumnView(1);
  EXPECT_DOUBLE_EQ(cc(0, 0), 10.0);
  EXPECT_DOUBLE_EQ(cc(1, 2), 15.0);
}

TEST(Array3D, RowVariable)
{
  tuvx::Array3D<double> a(2, 3, 1);
  auto tmp = a.GetRowVariable();

  tmp(0, 0) = 1.0;
  tmp(0, 1) = 2.0;
  tmp(0, 2) = 3.0;
  tmp(1, 0) = 4.0;
  tmp(1, 1) = 5.0;
  tmp(1, 2) = 6.0;

  EXPECT_DOUBLE_EQ(tmp(0, 0), 1.0);
  EXPECT_DOUBLE_EQ(tmp(1, 2), 6.0);
}

// ForEachRow: col 1 = col 0 * 2 across all (i, j) positions.
TEST(Array3D, ForEachRowScale)
{
  tuvx::Array3D<double> a(2, 3, 2);
  for (std::size_t i = 0; i < 2; ++i)
  {
    for (std::size_t j = 0; j < 3; ++j)
    {
      a(i, j, 0) = static_cast<double>((i * 3) + j + 1);
    }
  }

  a.ForEachRow(
      [](const double &src, double &dst) { dst = src * 2.0; },
      a.GetConstColumnView(0),
      a.GetColumnView(1));

  for (std::size_t i = 0; i < 2; ++i)
  {
    for (std::size_t j = 0; j < 3; ++j)
    {
      EXPECT_DOUBLE_EQ(a(i, j, 1), a(i, j, 0) * 2.0);
    }
  }
}

// ForEachRow with RowVariable as intermediate.
TEST(Array3D, ForEachRowWithRowVariable)
{
  tuvx::Array3D<double> a(2, 2, 3);  // [dim1=2, dim2=2, col=3]
  for (std::size_t i = 0; i < 2; ++i)
  {
    for (std::size_t j = 0; j < 2; ++j)
    {
      a(i, j, 0) = static_cast<double>((i * 2) + j + 1);    // 1,2,3,4
      a(i, j, 1) = static_cast<double>(((i * 2) + j) * 10);  // 0,10,20,30
    }
  }

  auto tmp = a.GetRowVariable();
  a.ForEachRow(
      [](const double &x, const double &y, double &t) { t = x + y; },
      a.GetConstColumnView(0),
      a.GetConstColumnView(1),
      tmp);

  a.ForEachRow(
      [](const double &t, double &out) { out = t * 2.0; },
      tmp,
      a.GetColumnView(2));

  EXPECT_DOUBLE_EQ(a(0, 0, 2), 2.0);   // (1 + 0) * 2
  EXPECT_DOUBLE_EQ(a(0, 1, 2), 24.0);  // (2 + 10) * 2
  EXPECT_DOUBLE_EQ(a(1, 0, 2), 46.0);  // (3 + 20) * 2
  EXPECT_DOUBLE_EQ(a(1, 1, 2), 68.0);  // (4 + 30) * 2
}

// Function: build a reusable operation and apply it to two arrays.
TEST(Array3D, Function)
{
  tuvx::Array3D<double> proto(2, 2, 3);

  auto func = tuvx::Array3D<double>::Function(
      [](auto &arr) {
        arr.ForEachRow(
            [](const double &src, double &dst) { dst = src * 5.0; },
            arr.GetConstColumnView(0),
            arr.GetColumnView(1));
      },
      proto);

  tuvx::Array3D<double> a(2, 2, 3);
  a(0, 0, 0) = 1.0;
  a(0, 1, 0) = 2.0;
  a(1, 0, 0) = 3.0;
  a(1, 1, 0) = 4.0;
  func(a);
  EXPECT_DOUBLE_EQ(a(0, 0, 1), 5.0);
  EXPECT_DOUBLE_EQ(a(0, 1, 1), 10.0);
  EXPECT_DOUBLE_EQ(a(1, 0, 1), 15.0);
  EXPECT_DOUBLE_EQ(a(1, 1, 1), 20.0);

  tuvx::Array3D<double> b(2, 2, 3);
  b(0, 0, 0) = 10.0;
  b(0, 1, 0) = 20.0;
  b(1, 0, 0) = 30.0;
  b(1, 1, 0) = 40.0;
  func(b);
  EXPECT_DOUBLE_EQ(b(0, 0, 1), 50.0);
  EXPECT_DOUBLE_EQ(b(1, 1, 1), 200.0);
}

// ForEachRow with Array2D as view (broadcasts operator()(i, j) naturally).
TEST(Array3D, ForEachRowWithArray2D)
{
  // weights [wl=2, z=3, col=2], temperature [wl=2, z=3]
  // (temperature is constant across columns in this test — treat it as a 2D table)
  tuvx::Array3D<double> weights(2, 3, 2);
  tuvx::Array2D<double> temp_table(2, 3);  // [wl × z] — no column dim

  for (std::size_t i = 0; i < 2; ++i)
  {
    for (std::size_t j = 0; j < 3; ++j)
    {
      weights(i, j, 0) = 1.0;
      temp_table(i, j) = static_cast<double>((i * 3) + j + 1);  // 1..6
    }
  }

  // col 1 = col 0 * temp_table(i, j)
  weights.ForEachRow(
      [](const double &w, const double &t, double &out) { out = w * t; },
      weights.GetConstColumnView(0),
      temp_table,
      weights.GetColumnView(1));

  for (std::size_t i = 0; i < 2; ++i)
  {
    for (std::size_t j = 0; j < 3; ++j)
    {
      EXPECT_DOUBLE_EQ(weights(i, j, 1), temp_table(i, j));
    }
  }
}
