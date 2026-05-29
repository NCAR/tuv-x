// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#include <tuvx/util/array2d.hpp>

#include <tuvx/util/array1d.hpp>

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

TEST(Array2D, ValueType)
{
  static_assert(std::is_same_v<tuvx::Array2D<double>::value_type, double>);
  static_assert(std::is_same_v<tuvx::Array2D<int>::value_type, int>);
}

// [rows=3, cols=2]:
//   col 0: rows 0,1,2  -> values 1,3,5
//   col 1: rows 0,1,2  -> values 2,4,6
TEST(Array2D, ColumnView)
{
  tuvx::Array2D<double> a(3, 2);
  a(0, 0) = 1.0;
  a(0, 1) = 2.0;
  a(1, 0) = 3.0;
  a(1, 1) = 4.0;
  a(2, 0) = 5.0;
  a(2, 1) = 6.0;

  auto col0 = a.GetColumnView(0);
  EXPECT_DOUBLE_EQ(col0[0], 1.0);
  EXPECT_DOUBLE_EQ(col0[1], 3.0);
  EXPECT_DOUBLE_EQ(col0[2], 5.0);

  auto col1 = a.GetColumnView(1);
  EXPECT_DOUBLE_EQ(col1[0], 2.0);
  EXPECT_DOUBLE_EQ(col1[1], 4.0);
  EXPECT_DOUBLE_EQ(col1[2], 6.0);

  // Mutate through view
  col0[1] = 99.0;
  EXPECT_DOUBLE_EQ(a(1, 0), 99.0);
}

TEST(Array2D, ConstColumnView)
{
  tuvx::Array2D<double> a(3, 2);
  a(0, 0) = 1.0;
  a(1, 0) = 2.0;
  a(2, 0) = 3.0;

  const auto &ca = a;
  auto cv = ca.GetConstColumnView(0);
  EXPECT_DOUBLE_EQ(cv[0], 1.0);
  EXPECT_DOUBLE_EQ(cv[1], 2.0);
  EXPECT_DOUBLE_EQ(cv[2], 3.0);
}

// RowVariable stores per-row temporaries, same interface as a column view.
TEST(Array2D, RowVariable)
{
  tuvx::Array2D<double> a(3, 2);
  auto tmp = a.GetRowVariable();

  tmp[0] = 10.0;
  tmp[1] = 20.0;
  tmp[2] = 30.0;

  EXPECT_DOUBLE_EQ(tmp[0], 10.0);
  EXPECT_DOUBLE_EQ(tmp[1], 20.0);
  EXPECT_DOUBLE_EQ(tmp[2], 30.0);
}

// ForEachRow: copy column 0 into column 1, scaling by 2.
TEST(Array2D, ForEachRowScale)
{
  tuvx::Array2D<double> a(4, 2);
  a(0, 0) = 1.0;
  a(1, 0) = 2.0;
  a(2, 0) = 3.0;
  a(3, 0) = 4.0;

  a.ForEachRow(
      [](const double &src, double &dst) { dst = src * 2.0; },
      a.GetConstColumnView(0),
      a.GetColumnView(1));

  EXPECT_DOUBLE_EQ(a(0, 1), 2.0);
  EXPECT_DOUBLE_EQ(a(1, 1), 4.0);
  EXPECT_DOUBLE_EQ(a(2, 1), 6.0);
  EXPECT_DOUBLE_EQ(a(3, 1), 8.0);
  // Source column unchanged
  EXPECT_DOUBLE_EQ(a(0, 0), 1.0);
  EXPECT_DOUBLE_EQ(a(1, 0), 2.0);
}

// ForEachRow with RowVariable as intermediate storage.
TEST(Array2D, ForEachRowWithRowVariable)
{
  tuvx::Array2D<double> a(3, 3);
  a(0, 0) = 1.0;
  a(1, 0) = 2.0;
  a(2, 0) = 3.0;
  a(0, 1) = 10.0;
  a(1, 1) = 20.0;
  a(2, 1) = 30.0;

  // tmp = col0 + col1
  auto tmp = a.GetRowVariable();
  a.ForEachRow(
      [](const double &a0, const double &a1, double &t) { t = a0 + a1; },
      a.GetConstColumnView(0),
      a.GetConstColumnView(1),
      tmp);

  // col2 = tmp * 2
  a.ForEachRow(
      [](const double &t, double &out) { out = t * 2.0; },
      tmp,
      a.GetColumnView(2));

  EXPECT_DOUBLE_EQ(a(0, 2), 22.0);
  EXPECT_DOUBLE_EQ(a(1, 2), 44.0);
  EXPECT_DOUBLE_EQ(a(2, 2), 66.0);
}

// Function: build a reusable operation and apply it to two arrays.
TEST(Array2D, Function)
{
  tuvx::Array2D<double> proto(3, 2);

  // op: col1 = col0 * 3
  auto func = tuvx::Array2D<double>::Function(
      [](auto &arr) {
        arr.ForEachRow(
            [](const double &src, double &dst) { dst = src * 3.0; },
            arr.GetConstColumnView(0),
            arr.GetColumnView(1));
      },
      proto);

  tuvx::Array2D<double> a(3, 2);
  a(0, 0) = 1.0;
  a(1, 0) = 2.0;
  a(2, 0) = 3.0;
  func(a);
  EXPECT_DOUBLE_EQ(a(0, 1), 3.0);
  EXPECT_DOUBLE_EQ(a(1, 1), 6.0);
  EXPECT_DOUBLE_EQ(a(2, 1), 9.0);

  tuvx::Array2D<double> b(3, 2);
  b(0, 0) = 4.0;
  b(1, 0) = 5.0;
  b(2, 0) = 6.0;
  func(b);
  EXPECT_DOUBLE_EQ(b(0, 1), 12.0);
  EXPECT_DOUBLE_EQ(b(1, 1), 15.0);
  EXPECT_DOUBLE_EQ(b(2, 1), 18.0);
}

// ForEachRow with Array1D as a view (uses operator[](size_t)).
TEST(Array2D, ForEachRowWithArray1D)
{
  tuvx::Array2D<double> a(3, 2);

  tuvx::Array1D<double> offsets(3);
  offsets[0] = 100.0;
  offsets[1] = 200.0;
  offsets[2] = 300.0;

  a(0, 0) = 1.0;
  a(1, 0) = 2.0;
  a(2, 0) = 3.0;

  a.ForEachRow(
      [](const double &src, const double &offset, double &dst) { dst = src + offset; },
      a.GetConstColumnView(0),
      offsets,
      a.GetColumnView(1));

  EXPECT_DOUBLE_EQ(a(0, 1), 101.0);
  EXPECT_DOUBLE_EQ(a(1, 1), 202.0);
  EXPECT_DOUBLE_EQ(a(2, 1), 303.0);
}
