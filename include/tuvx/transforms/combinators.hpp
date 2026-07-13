// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
//
// Higher-order combinators that build composite TransformFuncs.
//
// Each combinator takes one or more TransformFunc values and returns a new
// TransformFunc that calculates composite weights.  They work with any
// TransformFunc - library factories, combinators, or user-written lambdas.
#pragma once

#include <tuvx/transforms/transform_func.hpp>

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <vector>

namespace tuvx
{

  // ---------------------------------------------------------------------------
  // in_region

  /// @brief Returns a TransformFunc that applies @p t only within [wl_min, wl_max].
  ///
  /// Weights are set to zero outside the specified wavelength range.
  ///
  /// @param wl_min  Lower bound (inclusive), meters.
  /// @param wl_max  Upper bound (inclusive), meters.
  /// @param t       Transform to apply within the region.
  template<typename ArrayPolicy = Array3D<double>>
  auto in_region(
      typename ArrayPolicy::value_type wl_min,
      typename ArrayPolicy::value_type wl_max,
      TransformFunc<ArrayPolicy> t) -> TransformFunc<ArrayPolicy>
  {
    return [wl_min, wl_max, t = std::move(t)](
               const AtmosphericState<ArrayPolicy>& state, Array3D<typename ArrayPolicy::value_type>& weights)
    {
      using T = typename ArrayPolicy::value_type;
      const auto n_wl = weights.Size1();
      const auto n_z = weights.Size2();
      const auto n_col = weights.Size3();

      Array3D<T> tmp(n_wl, n_z, n_col);
      t(state, tmp);

      for (std::size_t wl = 0; wl < n_wl; ++wl)
      {
        const auto lambda = state.wavelength_grid_.mid_points_(wl, 0);
        const T fill = (lambda >= wl_min && lambda <= wl_max) ? T{ 1 } : T{ 0 };
        for (std::size_t z = 0; z < n_z; ++z)
        {
          for (std::size_t col = 0; col < n_col; ++col)
          {
            weights(wl, z, col) = tmp(wl, z, col) * fill;
          }
        }
      }
    };
  }

  // ---------------------------------------------------------------------------
  // multiply

  /// @brief Returns a TransformFunc whose weights are @p a weights multiplied by @p b weights.
  ///
  /// Both transforms are evaluated into separate buffers; the results are
  /// combined element-wise.
  template<typename ArrayPolicy = Array3D<double>>
  auto multiply(TransformFunc<ArrayPolicy> a, TransformFunc<ArrayPolicy> b) -> TransformFunc<ArrayPolicy>
  {
    return [a = std::move(a), b = std::move(b)](
               const AtmosphericState<ArrayPolicy>& state, Array3D<typename ArrayPolicy::value_type>& weights)
    {
      using T = typename ArrayPolicy::value_type;
      const auto n_wl = weights.Size1();
      const auto n_z = weights.Size2();
      const auto n_col = weights.Size3();

      a(state, weights);
      Array3D<T> tmp(n_wl, n_z, n_col);
      b(state, tmp);

      auto it_w = weights.begin();
      for (const auto& v : tmp)
      {
        *it_w *= v;
        ++it_w;
      }
    };
  }

  // ---------------------------------------------------------------------------
  // add

  /// @brief Returns a TransformFunc whose weights are @p a weights plus @p b weights.
  template<typename ArrayPolicy = Array3D<double>>
  auto add(TransformFunc<ArrayPolicy> a, TransformFunc<ArrayPolicy> b) -> TransformFunc<ArrayPolicy>
  {
    return [a = std::move(a), b = std::move(b)](
               const AtmosphericState<ArrayPolicy>& state, Array3D<typename ArrayPolicy::value_type>& weights)
    {
      using T = typename ArrayPolicy::value_type;
      const auto n_wl = weights.Size1();
      const auto n_z = weights.Size2();
      const auto n_col = weights.Size3();

      a(state, weights);
      Array3D<T> tmp(n_wl, n_z, n_col);
      b(state, tmp);

      auto it_w = weights.begin();
      for (const auto& v : tmp)
      {
        *it_w += v;
        ++it_w;
      }
    };
  }

  // ---------------------------------------------------------------------------
  // piecewise

  /// @brief One region in a piecewise transform.
  template<typename ArrayPolicy = Array3D<double>>
  struct PiecewiseRegion
  {
    typename ArrayPolicy::value_type wl_min_{};
    typename ArrayPolicy::value_type wl_max_{};
    TransformFunc<ArrayPolicy> transform_{};
  };

  /// @brief Returns a TransformFunc that applies a different sub-transform per wavelength region.
  ///
  /// For each wavelength bin, the first region whose [wl_min, wl_max] contains
  /// the bin mid-point is applied; wavelengths not covered by any region get
  /// zero weight.
  ///
  /// @param regions  Ordered list of (wl_min, wl_max, TransformFunc) tuples.
  template<typename ArrayPolicy = Array3D<double>>
  auto piecewise(std::vector<PiecewiseRegion<ArrayPolicy>> regions) -> TransformFunc<ArrayPolicy>
  {
    return [regions = std::move(regions)](
               const AtmosphericState<ArrayPolicy>& state, Array3D<typename ArrayPolicy::value_type>& weights)
    {
      using T = typename ArrayPolicy::value_type;
      const auto n_wl = weights.Size1();
      const auto n_z = weights.Size2();
      const auto n_col = weights.Size3();

      // Zero all weights - regions cover only their range.
      for (auto& w : weights)
      {
        w = T{ 0 };
      }

      Array3D<T> tmp(n_wl, n_z, n_col);
      for (const auto& region : regions)
      {
        region.transform_(state, tmp);

        for (std::size_t wl = 0; wl < n_wl; ++wl)
        {
          const auto lambda = state.wavelength_grid_.mid_points_(wl, 0);
          if (lambda >= region.wl_min_ && lambda <= region.wl_max_)
          {
            for (std::size_t z = 0; z < n_z; ++z)
            {
              for (std::size_t col = 0; col < n_col; ++col)
              {
                weights(wl, z, col) = tmp(wl, z, col);
              }
            }
          }
        }
      }
    };
  }

  // ---------------------------------------------------------------------------
  // clamp

  /// @brief Returns a TransformFunc that clamps weights to [min_val, max_val].
  ///
  /// @param t        Source transform.
  /// @param min_val  Minimum weight value.
  /// @param max_val  Maximum weight value.
  template<typename ArrayPolicy = Array3D<double>>
  auto clamp(
      TransformFunc<ArrayPolicy> t,
      typename ArrayPolicy::value_type min_val,
      typename ArrayPolicy::value_type max_val) -> TransformFunc<ArrayPolicy>
  {
    return [t = std::move(t), min_val, max_val](
               const AtmosphericState<ArrayPolicy>& state, Array3D<typename ArrayPolicy::value_type>& weights)
    {
      t(state, weights);
      for (auto& w : weights)
      {
        w = std::clamp(w, min_val, max_val);
      }
    };
  }

  // ---------------------------------------------------------------------------
  // override_band
  //
  // Named spectral bands and their wavelength limits:
  //   Lyman-alpha:    [121.4, 121.9] nm
  //   Schumann-Runge: [175.4, 206.2] nm

  /// @brief A named wavelength band with known spectral limits (meters).
  struct SpectralBand
  {
    double wl_min_{};
    double wl_max_{};
  };

  /// @brief Look up a named spectral band by its conventional name.
  ///
  /// Supported names: "lyman-alpha", "schumann-runge".
  /// Throws std::invalid_argument for unknown names.
  inline SpectralBand named_band(std::string_view name)
  {
    if (name == "lyman-alpha")
    {
      return { .wl_min_ = 1.214e-7, .wl_max_ = 1.219e-7 };
    }
    if (name == "schumann-runge")
    {
      return { .wl_min_ = 1.754e-7, .wl_max_ = 2.062e-7 };
    }
    throw std::invalid_argument("Unknown spectral band: " + std::string(name));
  }

  /// @brief Returns a TransformFunc that replaces weights in a named band with @p value.
  ///
  /// The base transform @p t is evaluated first; then every weight whose
  /// wavelength bin mid-point falls within the named band is overwritten with
  /// @p value.
  ///
  /// @param name   Band name: "lyman-alpha" or "schumann-runge".
  /// @param value  Override weight value within the band.
  /// @param t      Base transform evaluated before the override.
  template<typename ArrayPolicy = Array3D<double>>
  auto override_band(std::string_view name, typename ArrayPolicy::value_type value, TransformFunc<ArrayPolicy> t)
      -> TransformFunc<ArrayPolicy>
  {
    const SpectralBand band = named_band(name);  // validate at factory time
    return [band, value, t = std::move(t)](
               const AtmosphericState<ArrayPolicy>& state, Array3D<typename ArrayPolicy::value_type>& weights)
    {
      t(state, weights);
      const auto n_wl = weights.Size1();
      const auto n_z = weights.Size2();
      const auto n_col = weights.Size3();
      for (std::size_t wl = 0; wl < n_wl; ++wl)
      {
        const auto lambda = state.wavelength_grid_.mid_points_(wl, 0);
        if (lambda >= band.wl_min_ && lambda <= band.wl_max_)
        {
          for (std::size_t z = 0; z < n_z; ++z)
          {
            for (std::size_t col = 0; col < n_col; ++col)
            {
              weights(wl, z, col) = value;
            }
          }
        }
      }
    };
  }

  // ---------------------------------------------------------------------------
  // scale

  /// @brief Returns a TransformFunc whose weights are @p factor times those of @p t.
  ///
  /// @param factor  Scalar multiplier.
  /// @param t       Source transform.
  template<typename ArrayPolicy = Array3D<double>>
  auto scale(typename ArrayPolicy::value_type factor, TransformFunc<ArrayPolicy> t) -> TransformFunc<ArrayPolicy>
  {
    return [factor, t = std::move(t)](
               const AtmosphericState<ArrayPolicy>& state, Array3D<typename ArrayPolicy::value_type>& weights)
    {
      t(state, weights);
      for (auto& w : weights)
      {
        w *= factor;
      }
    };
  }

  // ---------------------------------------------------------------------------
  // normalize

  /// @brief Returns a TransformFunc whose weights are @p t's weights divided by
  ///        their sum over the wavelength axis (per height/column).
  ///
  /// This is a reduction: it evaluates @p t, sums the weights over wavelength for
  /// each (height, column), and divides each weight by that sum. A column whose
  /// weights sum to zero is left unchanged (avoiding a divide by zero).
  ///
  /// Used for normalized action spectra / cross-sections (e.g. a Gaussian
  /// spectral weight normalized to unit sum over the model grid).
  template<typename ArrayPolicy = Array3D<double>>
  auto normalize(TransformFunc<ArrayPolicy> t) -> TransformFunc<ArrayPolicy>
  {
    return [t = std::move(t)](const AtmosphericState<ArrayPolicy>& state, Array3D<typename ArrayPolicy::value_type>& weights)
    {
      using T = typename ArrayPolicy::value_type;
      t(state, weights);
      const auto n_wl = weights.Size1();
      const auto n_z = weights.Size2();
      const auto n_col = weights.Size3();
      for (std::size_t z = 0; z < n_z; ++z)
      {
        for (std::size_t col = 0; col < n_col; ++col)
        {
          T sum = 0;
          for (std::size_t wl = 0; wl < n_wl; ++wl)
          {
            sum += weights(wl, z, col);
          }
          if (sum != T{ 0 })
          {
            for (std::size_t wl = 0; wl < n_wl; ++wl)
            {
              weights(wl, z, col) /= sum;
            }
          }
        }
      }
    };
  }

}  // namespace tuvx
