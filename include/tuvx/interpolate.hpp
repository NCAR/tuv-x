// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
//
// Grid interpolation primitives.
//
// These resample tabulated data given on a source axis onto a target axis; a
// typical use is matching a molecular cross-section, quantum yield, or action
// spectrum onto the model wavelength grid.  They are pure numerics - arrays in,
// arrays out, no I/O.
//
// Two axis conventions are in play:
//   * interpolate_linear treats both axes as discrete POINTS and returns one
//     value per target point (length == x_target size).
//   * the three area/overlap methods treat the axes as bin EDGES and return one
//     value per target BIN (length == x_target size - 1).  Conserving callers
//     therefore pass the grid edges, not the mid-points.
//
// None of the interpolators extrapolate.  To control behavior where the source
// does not span the target, pad the source with add_point before interpolating
// (fill with zero, the boundary value, or a constant).
#pragma once

#include <tuvx/util/array1d.hpp>

#include <algorithm>
#include <cstddef>
#include <functional>
#include <stdexcept>

namespace tuvx
{

  // ---------------------------------------------------------------------------
  // interpolate_linear

  /// @brief Standard point-to-point linear interpolation.
  ///
  /// Maps @p y_source, given on the discrete points @p x_source, onto the
  /// discrete target points @p x_target.  Target points that fall outside the
  /// source range are set to zero (no extrapolation).
  ///
  /// @param x_target  Target points (monotonically increasing).
  /// @param x_source  Source points (monotonically increasing).
  /// @param y_source  Source values, one per source point.
  /// @return Interpolated values, one per target point (length == x_target size).
  template<typename T = double>
  auto interpolate_linear(const Array1D<T>& x_target, const Array1D<T>& x_source, const Array1D<T>& y_source) -> Array1D<T>
  {
    const std::size_t n = x_source.Size();
    Array1D<T> y_target(x_target.Size());  // value-initialized to zero

    std::ranges::transform(
        x_target,
        y_target.begin(),
        [&](T x) -> T
        {
          if (n < 2 || x < x_source[0] || x >= x_source[n - 1])
          {
            return T{ 0 };
          }
          // upper_bound gives the first source point strictly greater than x, so
          // the bracketing segment is [j, j + 1] with x in [x_source[j], x_source[j + 1]).
          const std::size_t j = static_cast<std::size_t>(std::ranges::upper_bound(x_source, x) - x_source.begin()) - 1;
          const T slope = (y_source[j + 1] - y_source[j]) / (x_source[j + 1] - x_source[j]);
          return y_source[j] + (slope * (x - x_source[j]));
        });
    return y_target;
  }

  // ---------------------------------------------------------------------------
  // interpolate_conserving

  /// @brief Area-conserving linear interpolation onto target bins.
  ///
  /// Maps @p y_source, given on the discrete points @p x_source, onto the target
  /// bins defined by the edges @p x_target.  The value in each bin is the average
  /// of the piecewise-linear source curve over that bin (the trapezoidal area
  /// under the curve divided by the bin width).  This is the standard choice for
  /// re-gridding cross sections, quantum yields, and action spectra.
  ///
  /// The source must fully span the target range; extrapolation is not
  /// permitted.  Throws std::invalid_argument if either axis is not
  /// monotonically increasing, or if the source does not overlap the target.
  ///
  /// @param x_target  Target bin edges (monotonically increasing), length nbins+1.
  /// @param x_source  Source points (monotonically increasing).
  /// @param y_source  Source values, one per source point.
  /// @return Bin averages, one per target bin (length == x_target size - 1).
  template<typename T = double>
  auto interpolate_conserving(const Array1D<T>& x_target, const Array1D<T>& x_source, const Array1D<T>& y_source)
      -> Array1D<T>
  {
    const std::size_t n = x_source.Size();
    const std::size_t ng = x_target.Size();

    if (std::ranges::adjacent_find(x_source, std::greater_equal<T>{}) != x_source.end())
    {
      throw std::invalid_argument("interpolate_conserving: source grid must be monotonically increasing");
    }
    if (std::ranges::adjacent_find(x_target, std::greater_equal<T>{}) != x_target.end())
    {
      throw std::invalid_argument("interpolate_conserving: target grid must be monotonically increasing");
    }
    if (x_source[0] > x_target[0] || x_source[n - 1] < x_target[ng - 1])
    {
      throw std::invalid_argument("interpolate_conserving: source and target grids do not overlap");
    }

    Array1D<T> y_target(ng - 1);  // value-initialized to zero
    std::size_t start = 0;        // scan hint: first source segment that may reach this bin
    for (std::size_t i = 0; i + 1 < ng; ++i)
    {
      T area = 0;
      const T lo = x_target[i];
      const T hi = x_target[i + 1];
      std::size_t k = start;
      // skip source segments that lie entirely below this bin
      while (k + 1 < n && x_source[k + 1] <= lo)
      {
        start = (k == 0) ? 0 : k - 1;
        ++k;
      }
      // accumulate the trapezoidal area of each source segment's overlap with [lo, hi]
      while (k + 1 < n && x_source[k] < hi)
      {
        start = (k == 0) ? 0 : k - 1;
        const T a1 = std::max(x_source[k], lo);
        const T a2 = std::min(x_source[k + 1], hi);
        if (x_source[k + 1] != x_source[k])
        {
          const T slope = (y_source[k + 1] - y_source[k]) / (x_source[k + 1] - x_source[k]);
          const T b1 = y_source[k] + (slope * (a1 - x_source[k]));
          const T b2 = y_source[k] + (slope * (a2 - x_source[k]));
          area += T{ 0.5 } * (a2 - a1) * (b2 + b1);
        }
        ++k;
      }
      y_target[i] = area / (hi - lo);
    }
    return y_target;
  }

  // ---------------------------------------------------------------------------
  // interpolate_fractional_source

  /// @brief Re-bin binned data by fractional overlap, normalized to the source
  ///        bin width.
  ///
  /// The source data represent an integral quantity over each source bin (e.g.
  /// extra-terrestrial flux per wavelength interval).  Each target bin receives
  /// the sum of the fractional source-bin areas that overlap it, so the total
  /// integrated quantity is conserved.  Both axes are bin edges.
  ///
  /// Source bins outside the target range contribute zero.  If @p fold_in is
  /// true, source data extending past the last target edge are integrated and
  /// folded back into the last target bin (recommended when re-gridding profiles
  /// that set total optical depth).
  ///
  /// @param x_target  Target bin edges (monotonically increasing), length nbins+1.
  /// @param x_source  Source bin edges (monotonically increasing), length nsrc+1.
  /// @param y_source  Source values, one per source bin (length nsrc).
  /// @param fold_in   Fold "overhang" past the last target edge into the last bin.
  /// @return One value per target bin (length == x_target size - 1).
  template<typename T = double>
  auto interpolate_fractional_source(
      const Array1D<T>& x_target,
      const Array1D<T>& x_source,
      const Array1D<T>& y_source,
      bool fold_in = false) -> Array1D<T>
  {
    const std::size_t nto = x_target.Size();
    const std::size_t nfrom = x_source.Size();
    const std::size_t n_bins = nto - 1;
    Array1D<T> y_target(n_bins);  // value-initialized to zero

    std::size_t start = 0;  // scan hint: first source bin that may reach this target bin
    std::size_t j = 0;
    for (std::size_t i = 0; i + 1 < nto; ++i)
    {
      T sum = 0;
      j = start;
      // skip source bins entirely below this target bin
      while (j + 1 < nfrom && x_source[j + 1] < x_target[i])
      {
        start = j;
        ++j;
      }
      // add the fractional overlap of each source bin, normalized to source width
      while (j + 1 < nfrom && x_source[j] <= x_target[i + 1])
      {
        const T a1 = std::max(x_source[j], x_target[i]);
        const T a2 = std::min(x_source[j + 1], x_target[i + 1]);
        sum += (y_source[j] * (a2 - a1)) / (x_source[j + 1] - x_source[j]);
        ++j;
      }
      y_target[i] = sum;
    }

    // integrate the "overhang" past the last target edge and fold it into the last bin
    if (fold_in && j >= 1)
    {
      const T a1 = x_target[nto - 1];
      const T a2 = x_source[j];
      if (a2 > a1 || j + 1 < nfrom)
      {
        T tail = (y_source[j - 1] * (a2 - a1)) / (x_source[j] - x_source[j - 1]);
        for (std::size_t k = j; k + 1 < nfrom; ++k)
        {
          tail += y_source[k] * (x_source[k + 1] - x_source[k]);
        }
        y_target[n_bins - 1] += tail;
      }
    }
    return y_target;
  }

  // ---------------------------------------------------------------------------
  // interpolate_fractional_target

  /// @brief Re-bin binned data by fractional overlap, normalized to the target
  ///        bin width.
  ///
  /// Like interpolate_fractional_source, but each target bin's sum of overlapping
  /// source areas is divided by the target bin width, giving a bin average rather
  /// than a conserved integral.  Both axes are bin edges.
  ///
  /// @param x_target  Target bin edges (monotonically increasing), length nbins+1.
  /// @param x_source  Source bin edges (monotonically increasing), length nsrc+1.
  /// @param y_source  Source values, one per source bin (length nsrc).
  /// @param fold_in   Fold "overhang" past the last target edge into the last bin.
  /// @return One value per target bin (length == x_target size - 1).
  template<typename T = double>
  auto interpolate_fractional_target(
      const Array1D<T>& x_target,
      const Array1D<T>& x_source,
      const Array1D<T>& y_source,
      bool fold_in = false) -> Array1D<T>
  {
    const std::size_t n = x_source.Size();
    const std::size_t ng = x_target.Size();
    Array1D<T> y_target(ng - 1);  // value-initialized to zero

    std::size_t start = 0;  // scan hint: first source bin that may reach this target bin
    std::size_t j = 0;
    for (std::size_t i = 0; i + 1 < ng; ++i)
    {
      T sum = 0;
      j = start;
      // skip source bins entirely below this target bin
      while (j + 1 < n && x_source[j + 1] < x_target[i])
      {
        start = j;
        ++j;
      }
      // add the fractional overlap of each source bin; normalized below to target width
      while (j + 1 < n && x_source[j] <= x_target[i + 1])
      {
        const T a1 = std::max(x_source[j], x_target[i]);
        const T a2 = std::min(x_source[j + 1], x_target[i + 1]);
        sum += y_source[j] * (a2 - a1);
        ++j;
      }
      y_target[i] = sum / (x_target[i + 1] - x_target[i]);
    }

    // integrate the "overhang" past the last target edge and fold it into the last bin
    if (fold_in && j >= 1)
    {
      const T a1 = x_target[ng - 1];
      const T a2 = x_source[j];
      if (a2 > a1 || j + 1 < n)
      {
        T tail = (y_source[j - 1] * (a2 - a1)) / (x_source[j] - x_source[j - 1]);
        for (std::size_t k = j; k + 1 < n; ++k)
        {
          tail += y_source[k] * (x_source[k + 1] - x_source[k]);
        }
        y_target[ng - 2] += tail;
      }
    }
    return y_target;
  }

  // ---------------------------------------------------------------------------
  // add_point

  /// @brief Insert a single (xnew, ynew) point into ascending gridded data.
  ///
  /// @p x must be strictly ascending, and @p xnew must not equal an existing
  /// grid value.  The point is inserted so @p x stays ascending, and @p y is
  /// updated in lock-step.  Callers use this to pad a source table (e.g. with
  /// zeros or boundary values at x = 0 and x = 1e38) so it spans the target grid
  /// before an area-conserving interpolation, which is how extrapolation
  /// behavior is expressed.
  ///
  /// Throws std::invalid_argument if @p x is not strictly ascending or if
  /// @p xnew coincides with an existing grid value.
  ///
  /// @param x     Grid values, ascending (modified in place).
  /// @param y     Data values paired with @p x (modified in place).
  /// @param xnew  Grid value to insert.
  /// @param ynew  Data value to insert.
  template<typename T = double>
  void add_point(Array1D<T>& x, Array1D<T>& y, T xnew, T ynew)
  {
    auto& xv = x.AsVector();
    auto& yv = y.AsVector();

    if (std::ranges::adjacent_find(xv, std::greater_equal<T>{}) != xv.end())
    {
      throw std::invalid_argument("add_point: grid not monotonically increasing");
    }
    if (std::ranges::find(xv, xnew) != xv.end())
    {
      throw std::invalid_argument("add_point: xnew exactly matches an existing grid value");
    }

    const auto position = std::ranges::upper_bound(xv, xnew);
    const auto index = position - xv.begin();
    xv.insert(position, xnew);
    yv.insert(yv.begin() + index, ynew);
  }

}  // namespace tuvx
