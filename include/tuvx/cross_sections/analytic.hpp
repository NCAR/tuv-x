// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
//
// Purely analytic (temperature- and height-independent) cross-section
// factories.  Each returns a TransformFunc that broadcasts a wavelength-only
// formula over all height levels and columns.
//
// All output weights are in m² (SI).  Formulas are ported from:
//   main:src/cross_sections/rayliegh.F90
//   main:src/cross_sections/hobr-oh_br.F90
//   main:src/cross_sections/t_butyl_nitrate.F90
//   main:src/cross_sections/nitroxy_acetone.F90
//   main:src/cross_sections/nitroxy_ethanol.F90
#pragma once

#include <tuvx/transforms/combinators.hpp>
#include <tuvx/transforms/factories.hpp>

#include <cmath>

namespace tuvx::cross_sections
{

  // ---------------------------------------------------------------------------
  // rayleigh
  //
  // Nicolet (1984) Rayleigh scattering cross-section.
  //
  // Formula (from rayliegh.F90):
  //   mu  = lambda_nm / 1000   (micrometers)
  //   pwr = 3.6772 + 0.389*mu + 0.09426/mu   (mu <= 0.55)
  //       = 4.04                               (mu > 0.55)
  //   sigma = 4.02e-28 / mu^pwr   [cm²]
  //
  // No wavelength-range restriction; applies across the full model grid.

  /// @brief Returns a TransformFunc for Rayleigh scattering cross-section.
  ///
  /// The Nicolet (1984) parameterisation:
  /// \f[ \sigma(\lambda) = \frac{4.02 \times 10^{-28}}{\mu^{p(\mu)}} \f]
  /// where \f$\mu = \lambda\,[\text{nm}]/1000\f$ (micrometers) and
  /// \f$p(\mu) = 3.6772 + 0.389\mu + 0.09426/\mu\f$ for
  /// \f$\mu \le 0.55\f$, else \f$p = 4.04\f$.
  ///
  /// @tparam ArrayPolicy  Storage policy (default: Array3D<double>).
  template<typename ArrayPolicy = Array3D<double>>
  auto rayleigh() -> TransformFunc<ArrayPolicy>
  {
    return tuvx::wrap_analytic<ArrayPolicy>(
        [](typename ArrayPolicy::value_type lambda_m) -> typename ArrayPolicy::value_type
        {
          const auto mu = lambda_m * 1.0e6;  // m -> micrometers (Fortran: lambda_nm/1000)
          const auto pwr = (mu <= 0.55)
                               ? 3.6772 + 0.389 * mu + 0.09426 / mu
                               : 4.04;
          const auto sigma_cm2 = 4.02e-28 / std::pow(mu, pwr);
          return sigma_cm2 * 1.0e-4;  // cm² -> m²
        });
  }

  // ---------------------------------------------------------------------------
  // hobr
  //
  // Hypobromous acid (HOBr) cross-section.
  //
  // Formula (from hobr-oh_br.F90):
  //   sigma = [24.77 * exp(-109.80 * (ln(284.01/wl))^2)
  //          + 12.22 * exp( -93.63 * (ln(350.57/wl))^2)
  //          + 2.283 * exp(-242.40 * (ln(457.38/wl))^2)] * 1e-20  [cm²]
  //   where wl is wavelength in nm.
  //   Active for 250 <= wl_nm <= 550; zero elsewhere.

  /// @brief Returns a TransformFunc for HOBr cross-section.
  ///
  /// Sum of three log-normal functions (Burkholder et al.):
  /// \f[ \sigma(\lambda) = \left[
  ///     24.77 e^{-109.80 (\ln(284.01/\lambda_\text{nm}))^2}
  ///   + 12.22 e^{-93.63  (\ln(350.57/\lambda_\text{nm}))^2}
  ///   + 2.283 e^{-242.40 (\ln(457.38/\lambda_\text{nm}))^2}
  ///   \right] \times 10^{-20}\,\text{cm}^2 \f]
  /// Active in [250, 550] nm; zero outside.
  ///
  /// @tparam ArrayPolicy  Storage policy (default: Array3D<double>).
  template<typename ArrayPolicy = Array3D<double>>
  auto hobr() -> TransformFunc<ArrayPolicy>
  {
    return tuvx::in_region<ArrayPolicy>(
        250.0e-9,
        550.0e-9,
        tuvx::wrap_analytic<ArrayPolicy>(
            [](typename ArrayPolicy::value_type lambda_m) -> typename ArrayPolicy::value_type
            {
              const auto wl = lambda_m * 1.0e9;  // m -> nm
              const auto sigma_cm2 =
                  (24.77 * std::exp(-109.80 * std::pow(std::log(284.01 / wl), 2)) +
                   12.22 * std::exp(-93.63 * std::pow(std::log(350.57 / wl), 2)) +
                   2.283 * std::exp(-242.40 * std::pow(std::log(457.38 / wl), 2))) *
                  1.0e-20;
              return sigma_cm2 * 1.0e-4;  // cm² -> m²
            }));
  }

  // ---------------------------------------------------------------------------
  // nitrate_formula (internal helper)
  //
  // Shared exponential-quadratic formula used by t-butyl nitrate,
  // nitroxy acetone, and nitroxy ethanol (all from the Fortran):
  //
  //   sigma_cm2 = exp(c + wl_nm * (b + a * wl_nm))
  //   Active for wl_min_nm <= wl_nm <= wl_max_nm; zero elsewhere.

  namespace detail
  {
    template<typename ArrayPolicy>
    auto nitrate_formula(
        typename ArrayPolicy::value_type a,
        typename ArrayPolicy::value_type b,
        typename ArrayPolicy::value_type c,
        typename ArrayPolicy::value_type wl_min_m,
        typename ArrayPolicy::value_type wl_max_m) -> TransformFunc<ArrayPolicy>
    {
      return tuvx::in_region<ArrayPolicy>(
          wl_min_m,
          wl_max_m,
          tuvx::wrap_analytic<ArrayPolicy>(
              [a, b, c](typename ArrayPolicy::value_type lambda_m) -> typename ArrayPolicy::value_type
              {
                const auto wl = lambda_m * 1.0e9;  // m -> nm
                const auto sigma_cm2 = std::exp(c + wl * (b + a * wl));
                return sigma_cm2 * 1.0e-4;  // cm² -> m²
              }));
    }
  }  // namespace detail

  // ---------------------------------------------------------------------------
  // t_butyl_nitrate
  //
  // Constants from t_butyl_nitrate.F90:
  //   a = -0.993e-3, b = 0.5307, c = -115.5
  //   Active: 270 <= wl_nm <= 330

  /// @brief Returns a TransformFunc for tert-butyl nitrate cross-section.
  ///
  /// \f[ \sigma(\lambda) = \exp\!\left( c + \lambda_\text{nm}(b + a\,\lambda_\text{nm}) \right) \f]
  /// with \f$a = -0.993 \times 10^{-3}\f$, \f$b = 0.5307\f$, \f$c = -115.5\f$.
  /// Active in [270, 330] nm; zero outside.
  template<typename ArrayPolicy = Array3D<double>>
  auto t_butyl_nitrate() -> TransformFunc<ArrayPolicy>
  {
    return detail::nitrate_formula<ArrayPolicy>(
        -0.993e-3, 0.5307, -115.5, 270.0e-9, 330.0e-9);
  }

  // ---------------------------------------------------------------------------
  // nitroxy_acetone
  //
  // Constants from nitroxy_acetone.F90:
  //   a = -1.365e-3, b = 0.7834, c = -156.8
  //   Active: 284 <= wl_nm <= 335

  /// @brief Returns a TransformFunc for nitroxy acetone cross-section.
  ///
  /// Same formula as t_butyl_nitrate() with \f$a = -1.365 \times 10^{-3}\f$,
  /// \f$b = 0.7834\f$, \f$c = -156.8\f$.  Active in [284, 335] nm.
  template<typename ArrayPolicy = Array3D<double>>
  auto nitroxy_acetone() -> TransformFunc<ArrayPolicy>
  {
    return detail::nitrate_formula<ArrayPolicy>(
        -1.365e-3, 0.7834, -156.8, 284.0e-9, 335.0e-9);
  }

  // ---------------------------------------------------------------------------
  // nitroxy_ethanol
  //
  // Constants from nitroxy_ethanol.F90:
  //   a = -2.359e-3, b = 1.2478, c = -210.4
  //   Active: 270 <= wl_nm <= 306

  /// @brief Returns a TransformFunc for nitroxy ethanol cross-section.
  ///
  /// Same formula as t_butyl_nitrate() with \f$a = -2.359 \times 10^{-3}\f$,
  /// \f$b = 1.2478\f$, \f$c = -210.4\f$.  Active in [270, 306] nm.
  template<typename ArrayPolicy = Array3D<double>>
  auto nitroxy_ethanol() -> TransformFunc<ArrayPolicy>
  {
    return detail::nitrate_formula<ArrayPolicy>(
        -2.359e-3, 1.2478, -210.4, 270.0e-9, 306.0e-9);
  }

}  // namespace tuvx::cross_sections
