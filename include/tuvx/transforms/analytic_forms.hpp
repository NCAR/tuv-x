// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
//
// General, configurable analytic wavelength-only transform forms.
//
// These are species-agnostic mathematical shapes that recur across many
// cross-sections, quantum yields, and spectral weights.  All parameters are
// supplied by the caller through a named-field parameters struct (so every
// value carries a name at the call site); no constant is hard-coded and no
// species is named here.  Example configurations that bind specific parameter
// values live in tuvx/fixed_configuration.hpp.
//
// Each form evaluates its formula only within an optional [wl_min, wl_max]
// wavelength window and writes zero outside it.  Evaluating only inside the
// window matters for forms that can overflow (e.g. exp of a polynomial):
// composing a plain transform with the in_region() combinator would evaluate
// the inner formula everywhere and could produce inf*0 = NaN outside the band.
#pragma once

#include <tuvx/transforms/transform_func.hpp>

#include <cmath>
#include <cstddef>
#include <limits>
#include <utility>
#include <vector>

namespace tuvx
{

  namespace detail
  {
    // Build a TransformFunc from a per-wavelength function f(lambda_m).
    // f is evaluated only where wl_min <= lambda <= wl_max; the result is
    // broadcast across all height levels and columns.  Outside the window the
    // weight is zero and f is never called.
    template<typename ArrayPolicy, typename F>
    auto wavelength_form(F f, typename ArrayPolicy::value_type wl_min, typename ArrayPolicy::value_type wl_max)
        -> TransformFunc<ArrayPolicy>
    {
      return [f = std::move(f), wl_min, wl_max](
                 const AtmosphericState<ArrayPolicy> &state, Array3D<typename ArrayPolicy::value_type> &weights)
      {
        using T = typename ArrayPolicy::value_type;
        const auto n_wl = weights.Size1();
        const auto n_z = weights.Size2();
        const auto n_col = weights.Size3();
        for (std::size_t wl = 0; wl < n_wl; ++wl)
        {
          const auto lambda = state.wavelength_grid_.mid_points_(wl, 0);
          const T value = (lambda >= wl_min && lambda <= wl_max) ? static_cast<T>(f(lambda)) : T{ 0 };
          for (std::size_t z = 0; z < n_z; ++z)
          {
            for (std::size_t col = 0; col < n_col; ++col)
            {
              weights(wl, z, col) = value;
            }
          }
        }
      };
    }

    // Build a TransformFunc from a per-element function f(lambda_m, temperature_K).
    // Like wavelength_form, but f also receives the local temperature, so it is
    // evaluated per (wavelength, height, column) rather than broadcast.  f is
    // called only where wl_min <= lambda <= wl_max; outside the window the weight
    // is zero and f is never called (so a formula that would overflow outside its
    // band never gets the chance to).
    template<typename ArrayPolicy, typename F>
    auto wavelength_temperature_form(
        F f,
        typename ArrayPolicy::value_type wl_min,
        typename ArrayPolicy::value_type wl_max) -> TransformFunc<ArrayPolicy>
    {
      return [f = std::move(f), wl_min, wl_max](
                 const AtmosphericState<ArrayPolicy> &state,
                 Array3D<typename ArrayPolicy::value_type> &weights)
      {
        using T = typename ArrayPolicy::value_type;
        const auto n_wl = weights.Size1();
        const auto n_z = weights.Size2();
        const auto n_col = weights.Size3();
        for (std::size_t wl = 0; wl < n_wl; ++wl)
        {
          const auto lambda = state.wavelength_grid_.mid_points_(wl, 0);
          const bool in_band = (lambda >= wl_min && lambda <= wl_max);
          for (std::size_t z = 0; z < n_z; ++z)
          {
            for (std::size_t col = 0; col < n_col; ++col)
            {
              weights(wl, z, col) =
                  in_band ? static_cast<T>(f(lambda, state.temperature_(z, col))) : T{ 0 };
            }
          }
        }
      };
    }
  }  // namespace detail

  // ---------------------------------------------------------------------------
  // log_normal_bands

  /// @brief One log-normal absorption band.
  ///
  /// Contributes \f$ A\,\exp(-w\,\ln^2(c/\lambda)) \f$ to the weight, i.e. a
  /// Gaussian in \f$\ln\lambda\f$ centred at wavelength \f$c\f$.
  template<typename ArrayPolicy = Array3D<double>>
  struct LogNormalBand
  {
    /// Peak amplitude at the band centre (same units as the output weight).
    typename ArrayPolicy::value_type amplitude_{};
    /// Width coefficient on \f$\ln^2(c/\lambda)\f$ (dimensionless; larger = narrower).
    typename ArrayPolicy::value_type width_{};
    /// Band centre wavelength (m).
    typename ArrayPolicy::value_type center_{};
  };

  /// @brief Parameters for log_normal_bands().
  template<typename ArrayPolicy = Array3D<double>>
  struct LogNormalBandsParameters
  {
    /// Bands to sum.
    std::vector<LogNormalBand<ArrayPolicy>> bands_{};
    /// Lower wavelength bound, inclusive (m). Default: no lower bound.
    typename ArrayPolicy::value_type wl_min_ = std::numeric_limits<typename ArrayPolicy::value_type>::lowest();
    /// Upper wavelength bound, inclusive (m). Default: no upper bound.
    typename ArrayPolicy::value_type wl_max_ = std::numeric_limits<typename ArrayPolicy::value_type>::max();
  };

  /// @brief Returns a TransformFunc equal to a sum of log-normal bands.
  ///
  /// \f[ w(\lambda) = \sum_i A_i \exp\!\left(-w_i\,\ln^2(c_i/\lambda)\right) \f]
  /// evaluated for \f$ \lambda_\text{min} \le \lambda \le \lambda_\text{max} \f$
  /// and zero outside.
  ///
  /// @tparam ArrayPolicy  Storage policy (default: Array3D<double>).
  /// @param params  Bands and optional wavelength window.
  template<typename ArrayPolicy = Array3D<double>>
  auto log_normal_bands(LogNormalBandsParameters<ArrayPolicy> params) -> TransformFunc<ArrayPolicy>
  {
    return detail::wavelength_form<ArrayPolicy>(
        [bands = std::move(params.bands_)](typename ArrayPolicy::value_type lambda) -> typename ArrayPolicy::value_type
        {
          using T = typename ArrayPolicy::value_type;
          T sum = 0;
          for (const auto &band : bands)
          {
            const T l = std::log(band.center_ / lambda);
            sum += band.amplitude_ * std::exp(-band.width_ * l * l);
          }
          return sum;
        },
        params.wl_min_,
        params.wl_max_);
  }

  // ---------------------------------------------------------------------------
  // exp_polynomial

  /// @brief Parameters for exp_polynomial().
  ///
  /// Computes \f$ s_\text{out}\,\exp\!\left(\sum_n c_n (s_\lambda \lambda)^n\right) \f$.
  template<typename ArrayPolicy = Array3D<double>>
  struct ExpPolynomialParameters
  {
    /// Polynomial coefficients, lowest order first (c0, c1, c2, ...).
    std::vector<typename ArrayPolicy::value_type> coefficients_{};
    /// Factor applied to lambda (m) before evaluating the polynomial
    /// (e.g. 1e9 to evaluate in nm while the grid is in m).
    typename ArrayPolicy::value_type wavelength_scale_ = 1;
    /// Factor applied to the exponential result (e.g. 1e-4 for cm^2 -> m^2).
    typename ArrayPolicy::value_type output_scale_ = 1;
    /// Lower wavelength bound, inclusive (m). Default: no lower bound.
    typename ArrayPolicy::value_type wl_min_ = std::numeric_limits<typename ArrayPolicy::value_type>::lowest();
    /// Upper wavelength bound, inclusive (m). Default: no upper bound.
    typename ArrayPolicy::value_type wl_max_ = std::numeric_limits<typename ArrayPolicy::value_type>::max();
  };

  /// @brief Returns a TransformFunc equal to the exponential of a wavelength polynomial.
  ///
  /// \f[ w(\lambda) = s_\text{out}\,
  ///     \exp\!\left( \sum_{n} c_n \, (s_\lambda\,\lambda)^n \right) \f]
  /// evaluated for \f$ \lambda_\text{min} \le \lambda \le \lambda_\text{max} \f$
  /// and zero outside.  The polynomial is evaluated by Horner's method.
  ///
  /// @tparam ArrayPolicy  Storage policy (default: Array3D<double>).
  /// @param params  Coefficients, scales, and optional wavelength window.
  template<typename ArrayPolicy = Array3D<double>>
  auto exp_polynomial(ExpPolynomialParameters<ArrayPolicy> params) -> TransformFunc<ArrayPolicy>
  {
    return detail::wavelength_form<ArrayPolicy>(
        [coeffs = std::move(params.coefficients_),
         wavelength_scale = params.wavelength_scale_,
         output_scale = params.output_scale_](typename ArrayPolicy::value_type lambda) -> typename ArrayPolicy::value_type
        {
          using T = typename ArrayPolicy::value_type;
          const T x = lambda * wavelength_scale;
          T poly = 0;
          // Horner: iterate from highest order down to c0.
          for (auto it = coeffs.rbegin(); it != coeffs.rend(); ++it)
          {
            poly = (poly * x) + *it;
          }
          return output_scale * std::exp(poly);
        },
        params.wl_min_,
        params.wl_max_);
  }

  // ---------------------------------------------------------------------------
  // bounded_analytic

  /// @brief Returns a TransformFunc that evaluates a user formula f(lambda_m, temperature_K).
  ///
  /// The result is \f$ f(\lambda, T) \f$ for
  /// \f$ \lambda_\text{min} \le \lambda \le \lambda_\text{max} \f$ and zero
  /// outside.  Unlike wrap_analytic(), the formula is evaluated only inside the
  /// wavelength window, so a formula that would overflow outside its active band
  /// (e.g. exp of a polynomial in wavelength) never gets evaluated there.
  ///
  /// Use this for temperature-dependent analytic cross-sections whose closed
  /// form is unique enough not to warrant a dedicated parameterised form; the
  /// formula and its constants live in the caller (e.g. fixed_configuration).
  ///
  /// @tparam ArrayPolicy  Storage policy (default: Array3D<double>).
  /// @param f       Callable (lambda_m, temperature_K) -> weight value.
  /// @param wl_min  Lower wavelength bound, inclusive (m). Default: no lower bound.
  /// @param wl_max  Upper wavelength bound, inclusive (m). Default: no upper bound.
  template<typename ArrayPolicy = Array3D<double>, typename F>
  auto bounded_analytic(
      F f,
      typename ArrayPolicy::value_type wl_min = std::numeric_limits<typename ArrayPolicy::value_type>::lowest(),
      typename ArrayPolicy::value_type wl_max = std::numeric_limits<typename ArrayPolicy::value_type>::max())
      -> TransformFunc<ArrayPolicy>
  {
    return detail::wavelength_temperature_form<ArrayPolicy>(std::move(f), wl_min, wl_max);
  }

}  // namespace tuvx
