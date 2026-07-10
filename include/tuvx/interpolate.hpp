// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
//
// Grid interpolation primitives.
//
// These are faithful C++ ports of the four interpolators in the Fortran
// tuvx_interpolate module (src/interpolate.F90) plus the add_point data-prep
// helper (src/util.F90).  They resample tabulated data given on a source axis
// onto a target axis; a typical use is matching a molecular cross-section,
// quantum yield, or action spectrum onto the model wavelength grid.
//
// They are pure numerics - arrays in, arrays out, no I/O.  The names are fully
// descriptive of the interpolation method, mirroring the Fortran type names
// (interpolator_linear_t, interpolator_conserving_t, ...).
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
#include <stdexcept>
#include <string>

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
    // 1-based read helpers so the body mirrors the Fortran line for line.
    auto xs = [&](int i) -> T { return x_source[static_cast<std::size_t>(i - 1)]; };
    auto ys = [&](int i) -> T { return y_source[static_cast<std::size_t>(i - 1)]; };

    const int n = static_cast<int>(x_source.Size());
    const int nt = static_cast<int>(x_target.Size());
    Array1D<T> y_target(static_cast<std::size_t>(nt));  // value-initialized to zero

    int jsave = 1;
    for (int i = 1; i <= nt; ++i)
    {
      const T xt = x_target[static_cast<std::size_t>(i - 1)];
      int j = jsave;
      while (true)
      {
        if (xt < xs(j) || xt >= xs(j + 1))
        {
          j = j + 1;
          if (j >= n)
          {
            break;
          }
        }
        else
        {
          const T slope = (ys(j + 1) - ys(j)) / (xs(j + 1) - xs(j));
          y_target[static_cast<std::size_t>(i - 1)] = ys(j) + (slope * (xt - xs(j)));
          jsave = j;
          break;
        }
      }
    }
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
    auto xs = [&](int i) -> T { return x_source[static_cast<std::size_t>(i - 1)]; };
    auto ys = [&](int i) -> T { return y_source[static_cast<std::size_t>(i - 1)]; };

    const int n = static_cast<int>(x_source.Size());
    const int ng = static_cast<int>(x_target.Size());

    // Grids must be monotonically increasing.
    for (int i = 1; i < n; ++i)
    {
      if (xs(i) >= xs(i + 1))
      {
        throw std::invalid_argument("interpolate_conserving: source grid must be monotonically increasing");
      }
    }
    for (int i = 1; i < ng; ++i)
    {
      if (x_target[static_cast<std::size_t>(i - 1)] >= x_target[static_cast<std::size_t>(i)])
      {
        throw std::invalid_argument("interpolate_conserving: target grid must be monotonically increasing");
      }
    }
    // Source must overlap the full target range.
    if (xs(1) > x_target[0] || xs(n) < x_target[static_cast<std::size_t>(ng - 1)])
    {
      throw std::invalid_argument("interpolate_conserving: source and target grids do not overlap");
    }

    Array1D<T> y_target(static_cast<std::size_t>(ng - 1));  // value-initialized to zero
    int jstart = 1;
    for (int i = 1; i <= ng - 1; ++i)
    {
      T area = 0;
      const T xgl = x_target[static_cast<std::size_t>(i - 1)];
      const T xgu = x_target[static_cast<std::size_t>(i)];
      // jstart is only a scan hint; keep it in range (a slightly-too-early start
      // rescans harmlessly and yields the same result). Both loops guard k < n,
      // so an out-of-range start simply produces a zero-area bin.
      int k = jstart < 1 ? 1 : jstart;
      // discard source points that lie entirely before this bin
      while (k < n && xs(k + 1) <= xgl)
      {
        jstart = k - 1;
        k = k + 1;
      }
      // accumulate trapezoidal area over each source segment's overlap with [xgl, xgu]
      while (k < n && xs(k) < xgu)
      {
        jstart = k - 1;
        const T a1 = std::max(xs(k), xgl);
        const T a2 = std::min(xs(k + 1), xgu);
        T darea = 0;
        if (xs(k + 1) != xs(k))
        {
          const T slope = (ys(k + 1) - ys(k)) / (xs(k + 1) - xs(k));
          const T b1 = ys(k) + (slope * (a1 - xs(k)));
          const T b2 = ys(k) + (slope * (a2 - xs(k)));
          darea = T{ 0.5 } * (a2 - a1) * (b2 + b1);
        }
        area = area + darea;
        k = k + 1;
      }
      y_target[static_cast<std::size_t>(i - 1)] = area / (xgu - xgl);
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
    auto xt = [&](int i) -> T { return x_target[static_cast<std::size_t>(i - 1)]; };
    auto xf = [&](int i) -> T { return x_source[static_cast<std::size_t>(i - 1)]; };
    auto yf = [&](int i) -> T { return y_source[static_cast<std::size_t>(i - 1)]; };

    const int nto = static_cast<int>(x_target.Size());
    const int nfrom = static_cast<int>(x_source.Size());
    const int ntobins = nto - 1;
    Array1D<T> y_target(static_cast<std::size_t>(ntobins));  // value-initialized to zero

    int jstart = 1;
    int j = jstart;
    for (int i = 1; i <= ntobins; ++i)
    {
      T sum = 0;
      j = jstart;
      // skip source bins entirely below this target bin
      while (j < nfrom && xf(j + 1) < xt(i))
      {
        jstart = j;
        j = j + 1;
      }
      // add the fractional overlap of each source bin, normalized to source width
      while (j < nfrom && xf(j) <= xt(i + 1))
      {
        const T a1 = std::max(xf(j), xt(i));
        const T a2 = std::min(xf(j + 1), xt(i + 1));
        sum = sum + ((yf(j) * (a2 - a1)) / (xf(j + 1) - xf(j)));
        j = j + 1;
      }
      y_target[static_cast<std::size_t>(i - 1)] = sum;
    }

    if (fold_in && j >= 2)
    {
      j = j - 1;
      const T a1 = xt(nto);
      const T a2 = xf(j + 1);
      if (a2 > a1 || j + 1 < nfrom)
      {
        T tail = (yf(j) * (a2 - a1)) / (xf(j + 1) - xf(j));
        for (int k = j + 1; k <= nfrom - 1; ++k)
        {
          tail = tail + (yf(k) * (xf(k + 1) - xf(k)));
        }
        y_target[static_cast<std::size_t>(ntobins - 1)] += tail;
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
    auto xt = [&](int i) -> T { return x_target[static_cast<std::size_t>(i - 1)]; };
    auto xs = [&](int i) -> T { return x_source[static_cast<std::size_t>(i - 1)]; };
    auto ys = [&](int i) -> T { return y_source[static_cast<std::size_t>(i - 1)]; };

    const int n = static_cast<int>(x_source.Size());
    const int ng = static_cast<int>(x_target.Size());
    Array1D<T> y_target(static_cast<std::size_t>(ng - 1));  // value-initialized to zero

    int jstart = 1;
    int j = jstart;
    for (int i = 1; i <= ng - 1; ++i)
    {
      T sum = 0;
      j = jstart;
      // skip source bins entirely below this target bin
      while (j < n && xs(j + 1) < xt(i))
      {
        jstart = j;
        j = j + 1;
      }
      // add the fractional overlap of each source bin, normalized below to target width
      while (j < n && xs(j) <= xt(i + 1))
      {
        const T a1 = std::max(xs(j), xt(i));
        const T a2 = std::min(xs(j + 1), xt(i + 1));
        sum = sum + (ys(j) * (a2 - a1));
        j = j + 1;
      }
      y_target[static_cast<std::size_t>(i - 1)] = sum / (xt(i + 1) - xt(i));
    }

    if (fold_in && j >= 2)
    {
      j = j - 1;
      const T a1 = xt(ng);
      const T a2 = xs(j + 1);
      if (a2 > a1 || j + 1 < n)
      {
        T tail = (ys(j) * (a2 - a1)) / (xs(j + 1) - xs(j));
        for (int k = j + 1; k <= n - 1; ++k)
        {
          tail = tail + (ys(k) * (xs(k + 1) - xs(k)));
        }
        y_target[static_cast<std::size_t>(ng - 2)] += tail;
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
    const std::size_t n = xv.size();

    for (std::size_t i = 1; i < n; ++i)
    {
      if (xv[i] <= xv[i - 1])
      {
        throw std::invalid_argument("add_point: grid not monotonically increasing");
      }
    }
    for (std::size_t i = 0; i < n; ++i)
    {
      if (xv[i] == xnew)
      {
        throw std::invalid_argument("add_point: xnew exactly matches an existing grid value");
      }
    }

    std::size_t insert = n;  // default: append at end (xnew > x.back())
    if (n == 0 || xnew < xv[0])
    {
      insert = 0;
    }
    else if (xnew < xv[n - 1])
    {
      for (std::size_t i = 1; i < n; ++i)
      {
        if (xv[i] > xnew)
        {
          insert = i;
          break;
        }
      }
    }

    xv.insert(xv.begin() + static_cast<std::ptrdiff_t>(insert), xnew);
    yv.insert(yv.begin() + static_cast<std::ptrdiff_t>(insert), ynew);
  }

}  // namespace tuvx
