// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
//
// Unit tests for the grid interpolation primitives in tuvx/interpolate.hpp.
//
// Each case uses inputs whose interpolated result is analytic (constant and
// linear sources average exactly; trapezoid areas are exact), so the expected
// values are computed by hand. Tolerance: relative 1e-12.
#include <tuvx/interpolate.hpp>
#include <tuvx/util/array1d.hpp>

#include <gtest/gtest.h>

#include <cstddef>
#include <initializer_list>
#include <stdexcept>
#include <vector>

namespace
{
  tuvx::Array1D<double> Make(std::initializer_list<double> vals)
  {
    tuvx::Array1D<double> a(vals.size());
    std::size_t i = 0;
    for (double v : vals)
    {
      a[i++] = v;
    }
    return a;
  }

  void ExpectClose(const tuvx::Array1D<double>& got, std::initializer_list<double> expected)
  {
    ASSERT_EQ(got.Size(), expected.size());
    std::size_t i = 0;
    for (double exp : expected)
    {
      const double g = got[i];
      if (exp == 0.0)
      {
        EXPECT_DOUBLE_EQ(g, 0.0) << "at index " << i;
      }
      else
      {
        EXPECT_LT(std::abs((g - exp) / exp), 1.0e-12) << "at index " << i << " got=" << g << " expected=" << exp;
      }
      ++i;
    }
  }
}  // namespace

TEST(Interpolate, Linear)
{
  // Below-range -> 0, exact node reproduction, interior points, above-range -> 0.
  const auto y = tuvx::interpolate_linear(
      Make({ 0.5, 1.0, 1.5, 2.5, 3.9, 4.5 }), Make({ 1.0, 2.0, 3.0, 4.0 }), Make({ 10.0, 20.0, 30.0, 40.0 }));
  ExpectClose(y, { 0.0, 10.0, 15.0, 25.0, 39.0, 0.0 });
}

TEST(Interpolate, ConservingConstantSource)
{
  // A constant source averages to the same constant in every bin.
  const auto y =
      tuvx::interpolate_conserving(Make({ 1.0, 2.0, 3.0, 4.0 }), Make({ 1.0, 2.0, 3.0, 4.0 }), Make({ 5.0, 5.0, 5.0, 5.0 }));
  ExpectClose(y, { 5.0, 5.0, 5.0 });
}

TEST(Interpolate, ConservingLinearSource)
{
  // Source y = x; bin average over [a,b] is the midpoint (a+b)/2.
  const auto y =
      tuvx::interpolate_conserving(Make({ 1.0, 2.0, 3.0, 4.0 }), Make({ 1.0, 2.0, 3.0, 4.0 }), Make({ 1.0, 2.0, 3.0, 4.0 }));
  ExpectClose(y, { 1.5, 2.5, 3.5 });
}

TEST(Interpolate, ConservingPartialBin)
{
  // Single target bin [1.5, 2.5] over source y = x averages to 2.0.
  const auto y =
      tuvx::interpolate_conserving(Make({ 1.5, 2.5 }), Make({ 1.0, 2.0, 3.0, 4.0 }), Make({ 1.0, 2.0, 3.0, 4.0 }));
  ExpectClose(y, { 2.0 });
}

TEST(Interpolate, ConservingRejectsNonMonotonicSource)
{
  EXPECT_THROW(
      tuvx::interpolate_conserving(Make({ 1.0, 2.0 }), Make({ 1.0, 1.0, 3.0 }), Make({ 1.0, 2.0, 3.0 })),
      std::invalid_argument);
}

TEST(Interpolate, ConservingRejectsNonOverlappingGrids)
{
  // Target extends beyond the source range -> no overlap -> throw.
  EXPECT_THROW(
      tuvx::interpolate_conserving(Make({ 1.0, 2.0, 5.0 }), Make({ 1.0, 2.0, 3.0 }), Make({ 1.0, 2.0, 3.0 })),
      std::invalid_argument);
}

TEST(Interpolate, FractionalSourceVsTarget)
{
  // Same binned input; the two methods differ only in width normalization.
  const auto src =
      tuvx::interpolate_fractional_source(Make({ 0.0, 1.5, 3.0 }), Make({ 0.0, 1.0, 2.0, 3.0 }), Make({ 10.0, 20.0, 30.0 }));
  ExpectClose(src, { 20.0, 40.0 });  // conserves total input mass 60

  const auto tgt =
      tuvx::interpolate_fractional_target(Make({ 0.0, 1.5, 3.0 }), Make({ 0.0, 1.0, 2.0, 3.0 }), Make({ 10.0, 20.0, 30.0 }));
  ExpectClose(tgt, { 13.333333333333334, 26.666666666666668 });
}

TEST(Interpolate, FractionalSourceFoldIn)
{
  const auto edges = Make({ 0.0, 1.0, 2.0, 3.0, 4.0 });
  const auto vals = Make({ 10.0, 20.0, 30.0, 40.0 });
  ExpectClose(tuvx::interpolate_fractional_source(Make({ 0.0, 2.0 }), edges, vals, false), { 30.0 });
  // fold_in adds the overhang bins [2,3]=30 and [3,4]=40 back into the last bin.
  ExpectClose(tuvx::interpolate_fractional_source(Make({ 0.0, 2.0 }), edges, vals, true), { 100.0 });
}

TEST(Interpolate, FractionalTargetFoldIn)
{
  const auto edges = Make({ 0.0, 1.0, 2.0, 3.0, 4.0 });
  const auto vals = Make({ 10.0, 20.0, 30.0, 40.0 });
  ExpectClose(tuvx::interpolate_fractional_target(Make({ 0.0, 2.0 }), edges, vals, false), { 15.0 });
  // The folded overhang (70) is added raw to the last bin, not width-normalized.
  ExpectClose(tuvx::interpolate_fractional_target(Make({ 0.0, 2.0 }), edges, vals, true), { 85.0 });
}

TEST(Interpolate, AddPointFront)
{
  auto x = Make({ 1.0, 2.0, 3.0 });
  auto y = Make({ 10.0, 20.0, 30.0 });
  tuvx::add_point(x, y, 0.5, 5.0);
  ExpectClose(x, { 0.5, 1.0, 2.0, 3.0 });
  ExpectClose(y, { 5.0, 10.0, 20.0, 30.0 });
}

TEST(Interpolate, AddPointMiddle)
{
  auto x = Make({ 1.0, 2.0, 3.0 });
  auto y = Make({ 10.0, 20.0, 30.0 });
  tuvx::add_point(x, y, 2.5, 25.0);
  ExpectClose(x, { 1.0, 2.0, 2.5, 3.0 });
  ExpectClose(y, { 10.0, 20.0, 25.0, 30.0 });
}

TEST(Interpolate, AddPointEnd)
{
  auto x = Make({ 1.0, 2.0, 3.0 });
  auto y = Make({ 10.0, 20.0, 30.0 });
  tuvx::add_point(x, y, 4.0, 40.0);
  ExpectClose(x, { 1.0, 2.0, 3.0, 4.0 });
  ExpectClose(y, { 10.0, 20.0, 30.0, 40.0 });
}

TEST(Interpolate, AddPointRejectsDuplicateAndNonMonotonic)
{
  auto x = Make({ 1.0, 2.0, 3.0 });
  auto y = Make({ 10.0, 20.0, 30.0 });
  EXPECT_THROW(tuvx::add_point(x, y, 2.0, 99.0), std::invalid_argument);

  auto xbad = Make({ 1.0, 3.0, 2.0 });
  auto ybad = Make({ 10.0, 20.0, 30.0 });
  EXPECT_THROW(tuvx::add_point(xbad, ybad, 4.0, 40.0), std::invalid_argument);
}
