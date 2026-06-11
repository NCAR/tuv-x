// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
//
// Example hard-coded configuration.
//
// This header demonstrates how specific species are assembled from the
// general, species-agnostic transform forms (tuvx/transforms/analytic_forms.hpp
// and the factory/combinator library).  It is the one place where particular
// parameter values, species identities, and literature citations live -- the
// role a host model or MUSICA's configuration layer would otherwise fill.
//
// The core transform library contains no species names or magic numbers; this
// file is a copyable example, and the regression tests use it to reproduce the
// reference values from the original implementation.
//
// All cross-sections return weights in m^2 (SI).
#pragma once

#include <tuvx/transforms/analytic_forms.hpp>
#include <tuvx/transforms/factories.hpp>

#include <cmath>

namespace tuvx::fixed_configuration
{

  /// @brief Rayleigh scattering cross-section (m^2).
  ///
  /// \f[ \sigma(\lambda) = \frac{4.02\times10^{-28}}{\mu^{p(\mu)}}\,\text{cm}^2,
  ///     \quad \mu = \lambda_\text{nm}/1000 \f]
  /// with \f$p(\mu) = 3.6772 + 0.389\,\mu + 0.09426/\mu\f$ for
  /// \f$\mu \le 0.55\f$, else \f$p = 4.04\f$.  Converted to m^2 with a factor
  /// of \f$10^{-4}\f$.  The reciprocal piecewise exponent has no general form,
  /// so it is written directly with wrap_analytic.
  ///
  /// \rst
  /// From WMO (1985), originally :cite:`Nicolet1984`.
  /// \endrst
  template<typename ArrayPolicy = Array3D<double>>
  auto rayleigh() -> TransformFunc<ArrayPolicy>
  {
    return tuvx::wrap_analytic<ArrayPolicy>(
        [](typename ArrayPolicy::value_type lambda_m) -> typename ArrayPolicy::value_type
        {
          const auto mu = lambda_m * 1.0e6;  // m -> micrometers (= lambda_nm/1000)
          const auto pwr = (mu <= 0.55) ? 3.6772 + (0.389 * mu) + (0.09426 / mu) : 4.04;
          const auto sigma_cm2 = 4.02e-28 / std::pow(mu, pwr);
          return sigma_cm2 * 1.0e-4;  // cm^2 -> m^2
        });
  }

  /// @brief HOBr cross-section (m^2): three log-normal bands, active 250-550 nm.
  ///
  /// Bands peak near 284, 351, and 457 nm.  Amplitudes are the published
  /// peak cross-sections (10^-20 cm^2) converted to m^2 (x10^-24 overall).
  ///
  /// \rst
  /// Cross-section data from :cite:`Ingham1998`.
  /// \endrst
  template<typename ArrayPolicy = Array3D<double>>
  auto hobr() -> TransformFunc<ArrayPolicy>
  {
    using T = typename ArrayPolicy::value_type;
    // amplitude = (coefficient * 1e-20 cm^2) * 1e-4 (cm^2 -> m^2) = coefficient * 1e-24
    return tuvx::log_normal_bands<ArrayPolicy>(
        { .bands_ = { LogNormalBand<ArrayPolicy>{
                          .amplitude_ = T{ 24.77e-24 }, .width_ = T{ 109.80 }, .center_ = T{ 284.01e-9 } },
                      LogNormalBand<ArrayPolicy>{
                          .amplitude_ = T{ 12.22e-24 }, .width_ = T{ 93.63 }, .center_ = T{ 350.57e-9 } },
                      LogNormalBand<ArrayPolicy>{
                          .amplitude_ = T{ 2.283e-24 }, .width_ = T{ 242.40 }, .center_ = T{ 457.38e-9 } } },
          .wl_min_ = T{ 250.0e-9 },
          .wl_max_ = T{ 550.0e-9 } });
  }

  // The three organic nitrates share the form
  //   sigma_cm2 = exp(c + b*lambda_nm + a*lambda_nm^2),  active in [wl_min, wl_max].
  // Expressed via exp_polynomial with the polynomial evaluated in nm
  // (wavelength_scale = 1e9) and the cm^2 -> m^2 conversion as output_scale (1e-4).
  //
  // The log-quadratic-in-wavelength form is characteristic of Roberts & Fajer
  // (1989); the specific coefficients below have not been verified against a
  // primary source, so no citation is asserted.

  /// @brief tert-butyl nitrate cross-section (m^2), active 270-330 nm.
  template<typename ArrayPolicy = Array3D<double>>
  auto t_butyl_nitrate() -> TransformFunc<ArrayPolicy>
  {
    using T = typename ArrayPolicy::value_type;
    return tuvx::exp_polynomial<ArrayPolicy>({ .coefficients_ = { T{ -115.5 }, T{ 0.5307 }, T{ -0.993e-3 } },
                                               .wavelength_scale_ = T{ 1.0e9 },
                                               .output_scale_ = T{ 1.0e-4 },
                                               .wl_min_ = T{ 270.0e-9 },
                                               .wl_max_ = T{ 330.0e-9 } });
  }

  /// @brief Nitroxy acetone cross-section (m^2), active 284-335 nm.
  template<typename ArrayPolicy = Array3D<double>>
  auto nitroxy_acetone() -> TransformFunc<ArrayPolicy>
  {
    using T = typename ArrayPolicy::value_type;
    return tuvx::exp_polynomial<ArrayPolicy>({ .coefficients_ = { T{ -156.8 }, T{ 0.7834 }, T{ -1.365e-3 } },
                                               .wavelength_scale_ = T{ 1.0e9 },
                                               .output_scale_ = T{ 1.0e-4 },
                                               .wl_min_ = T{ 284.0e-9 },
                                               .wl_max_ = T{ 335.0e-9 } });
  }

  /// @brief Nitroxy ethanol cross-section (m^2), active 270-306 nm.
  template<typename ArrayPolicy = Array3D<double>>
  auto nitroxy_ethanol() -> TransformFunc<ArrayPolicy>
  {
    using T = typename ArrayPolicy::value_type;
    return tuvx::exp_polynomial<ArrayPolicy>({ .coefficients_ = { T{ -210.4 }, T{ 1.2478 }, T{ -2.359e-3 } },
                                               .wavelength_scale_ = T{ 1.0e9 },
                                               .output_scale_ = T{ 1.0e-4 },
                                               .wl_min_ = T{ 270.0e-9 },
                                               .wl_max_ = T{ 306.0e-9 } });
  }

}  // namespace tuvx::fixed_configuration
