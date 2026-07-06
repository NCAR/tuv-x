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
#include <tuvx/transforms/combinators.hpp>
#include <tuvx/transforms/factories.hpp>
#include <tuvx/util/array1d.hpp>
#include <tuvx/util/array2d.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <initializer_list>
#include <utility>
#include <vector>

namespace tuvx::fixed_configuration
{

  namespace detail
  {
    /// Evaluate a polynomial by Horner's method. Coefficients are highest order
    /// first (e.g. {c_n, ..., c_1, c_0}).
    template<typename T>
    T horner(std::initializer_list<T> coefficients_high_to_low, T x)
    {
      T acc = 0;
      for (T c : coefficients_high_to_low)
      {
        acc = (acc * x) + c;
      }
      return acc;
    }

    /// A tabulated (wavelength_m, value) table.
    template<typename T>
    using TabulatedTable = std::vector<std::pair<T, T>>;

    /// Look up a wavelength (m) in a (wavelength_m, value) table, returning the
    /// matching value or zero if the wavelength is absent. This exact-match
    /// lookup is a deliberate placeholder for the conserving interpolator
    /// (Phases 4-5): it is faithful only when the model grid coincides with the
    /// tabulated grid, which is how the hybrid regression grids are pinned.
    template<typename T>
    T tabulated_lookup(const TabulatedTable<T> &table, T lambda_m)
    {
      for (const auto &[wavelength, value] : table)
      {
        if (std::abs(wavelength - lambda_m) <= T{ 1.0e-12 } * wavelength)
        {
          return value;
        }
      }
      return T{ 0 };
    }

    /// Build an Array1D from literal values, each multiplied by @p scale
    /// (e.g. 1e-4 to convert a cm^2 table to m^2).
    template<typename T>
    tuvx::Array1D<T> scaled_array(std::initializer_list<double> values, double scale)
    {
      tuvx::Array1D<T> out(values.size());
      std::size_t i = 0;
      for (double v : values)
      {
        out[i++] = static_cast<T>(v * scale);
      }
      return out;
    }

    /// Build an [n x 2] coefficient array for exponential_scaling with a zero
    /// zeroth-order term and the given per-wavelength first-order coefficients,
    /// so exponential_scaling yields exp(coefficient * (T - T_ref)).
    template<typename T>
    tuvx::Array2D<T> first_order_exponent(std::initializer_list<double> coefficients)
    {
      tuvx::Array2D<T> out(coefficients.size(), 2);
      std::size_t i = 0;
      for (double c : coefficients)
      {
        out(i, 0) = T{ 0 };
        out(i, 1) = static_cast<T>(c);
        ++i;
      }
      return out;
    }
  }  // namespace detail

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

  /// @brief N2O -> N2 + O(1D) cross-section (m^2), active 173-240 nm.
  ///
  /// \f[ \sigma(\lambda, T) = \exp\!\left( A(\lambda) + (T_\text{adj} - 300)\,e^{B(\lambda)} \right) \f]
  /// with \f$T_\text{adj} = \mathrm{clamp}(T, 194, 320)\f$, \f$A\f$ a degree-4
  /// and \f$B\f$ a degree-3 polynomial in wavelength (nm).
  template<typename ArrayPolicy = Array3D<double>>
  auto n2o() -> TransformFunc<ArrayPolicy>
  {
    using T = typename ArrayPolicy::value_type;
    return tuvx::bounded_analytic<ArrayPolicy>(
        [](T lambda_m, T temperature) -> T
        {
          const T wl = lambda_m * T{ 1.0e9 };  // nm
          const T t_adj = std::max(T{ 194.0 }, std::min(temperature, T{ 320.0 }));
          const T a = detail::horner<T>(
              { T{ 2.520672e-7 }, T{ -1.777846e-4 }, T{ 4.301146e-2 }, T{ -4.071805 }, T{ 68.21023 } }, wl);
          const T b_poly = detail::horner<T>({ T{ -1.881058e-5 }, T{ 1.111572e-2 }, T{ -2.116255 }, T{ 123.4014 } }, wl);
          const T b = (t_adj - T{ 300.0 }) * std::exp(b_poly);
          return std::exp(a + b) * T{ 1.0e-4 };  // cm^2 -> m^2
        },
        T{ 173.0e-9 },
        T{ 240.0e-9 });
  }

  /// @brief Cl2 -> Cl + Cl cross-section (m^2).
  ///
  /// \f[ \sigma(\lambda, T) = 10^{-20}\,\sqrt{\alpha}\,
  ///     \left[ 27.3\,e^{-99\alpha \ln^2(329.5/\lambda)}
  ///          + 0.932\,e^{-91.5\alpha \ln^2(406.5/\lambda)} \right]\,\text{cm}^2 \f]
  /// with \f$\alpha = \tanh(402.7/T)\f$ written as \f$(e^{2\beta}-1)/(e^{2\beta}+1)\f$,
  /// \f$\beta = 402.7/T\f$, and \f$\lambda\f$ in nm.
  template<typename ArrayPolicy = Array3D<double>>
  auto cl2() -> TransformFunc<ArrayPolicy>
  {
    using T = typename ArrayPolicy::value_type;
    return tuvx::bounded_analytic<ArrayPolicy>(
        [](T lambda_m, T temperature) -> T
        {
          const T wl = lambda_m * T{ 1.0e9 };  // nm
          const T beta = T{ 402.7 } / temperature;
          const T bb = std::exp(beta);
          const T bbsq = bb * bb;
          const T alpha = (bbsq - T{ 1.0 }) / (bbsq + T{ 1.0 });
          const T l1 = std::log(T{ 329.5 } / wl);
          const T l2 = std::log(T{ 406.5 } / wl);
          const T ex1 = T{ 27.3 } * std::exp(T{ -99.0 } * alpha * l1 * l1);
          const T ex2 = T{ 0.932 } * std::exp(T{ -91.5 } * alpha * l2 * l2);
          return T{ 1.0e-20 } * std::sqrt(alpha) * (ex1 + ex2) * T{ 1.0e-4 };  // cm^2 -> m^2
        });
  }

  /// @brief H2O2 -> OH + OH cross-section (m^2). Hybrid: analytic in 260-350 nm,
  ///        tabulated below.
  ///
  /// In band: \f$ \sigma = (\chi\,A(\lambda) + (1-\chi)\,B(\lambda))\,10^{-21}\,\text{cm}^2 \f$
  /// with \f$\chi = 1/(1 + e^{-1265/T})\f$, \f$T\f$ clamped to [200, 400], \f$A\f$
  /// degree-7 and \f$B\f$ degree-4 polynomials in nm. Below the band the
  /// tabulated base cross-section is used. The Fortran source constant
  /// \f$A_0 = 6.4761\times10^{4}\f$ is authoritative (the bundled .nc header text
  /// shows a typo'd exponent).
  template<typename ArrayPolicy = Array3D<double>>
  auto h2o2() -> TransformFunc<ArrayPolicy>
  {
    using T = typename ArrayPolicy::value_type;
    // Tabulated base below the analytic band, from cross_section.h2o2-oh_oh.nc (cm^2 -> m^2).
    const std::vector<std::pair<T, T>> tabulated = { { T{ 190.0e-9 }, T{ 6.72e-19 } * T{ 1.0e-4 } },
                                                     { T{ 195.0e-9 }, T{ 5.64e-19 } * T{ 1.0e-4 } },
                                                     { T{ 200.0e-9 }, T{ 4.75e-19 } * T{ 1.0e-4 } },
                                                     { T{ 205.0e-9 }, T{ 4.08e-19 } * T{ 1.0e-4 } },
                                                     { T{ 210.0e-9 }, T{ 3.57e-19 } * T{ 1.0e-4 } } };
    return [tabulated](const AtmosphericState<ArrayPolicy> &state, Array3D<T> &weights)
    {
      const auto n_wl = weights.Size1();
      const auto n_z = weights.Size2();
      const auto n_col = weights.Size3();
      for (std::size_t wl = 0; wl < n_wl; ++wl)
      {
        const T lambda_m = state.wavelength_grid_.mid_points_(wl, 0);
        const T wl_nm = lambda_m * T{ 1.0e9 };
        const bool in_band = (wl_nm >= T{ 260.0 } && wl_nm < T{ 350.0 });
        for (std::size_t z = 0; z < n_z; ++z)
        {
          for (std::size_t col = 0; col < n_col; ++col)
          {
            if (in_band)
            {
              const T t_adj = std::max(T{ 200.0 }, std::min(state.temperature_(z, col), T{ 400.0 }));
              const T chi = T{ 1.0 } / (T{ 1.0 } + std::exp(T{ -1265.0 } / t_adj));
              const T sum_a = detail::horner<T>(
                  { T{ 1.5534675e-13 },
                    T{ -2.652014e-10 },
                    T{ 1.6878206e-7 },
                    T{ -4.035101e-5 },
                    T{ -4.4589016e-3 },
                    T{ 4.535649 },
                    T{ -9.2170972e2 },
                    T{ 6.4761e4 } },
                  wl_nm);
              const T sum_b = detail::horner<T>(
                  { T{ -1.0924e-7 }, T{ -3.0493e-5 }, T{ 1.1522e-1 }, T{ -5.1351e1 }, T{ 6.8123e3 } }, wl_nm);
              weights(wl, z, col) = ((chi * sum_a) + ((T{ 1.0 } - chi) * sum_b)) * T{ 1.0e-21 } * T{ 1.0e-4 };
            }
            else
            {
              weights(wl, z, col) = detail::tabulated_lookup<T>(tabulated, lambda_m);
            }
          }
        }
      }
    };
  }

  /// @brief CHBr3 cross-section (m^2). Hybrid: analytic where 290 < lambda < 340 nm
  ///        AND 210 < T < 300 K, tabulated otherwise.
  ///
  /// In the analytic regime:
  /// \f[ \sigma = \exp\!\left( (C_0 - C_1\lambda)(T_0 - T) - (C_2 + C_3\lambda) \right)\,\text{cm}^2 \f]
  /// with \f$\lambda\f$ in nm.
  template<typename ArrayPolicy = Array3D<double>>
  auto chbr3() -> TransformFunc<ArrayPolicy>
  {
    using T = typename ArrayPolicy::value_type;
    // Tabulated base, from cross_section.chbr3.nc (cm^2 -> m^2).
    const std::vector<std::pair<T, T>> tabulated = {
      { T{ 190.0e-9 }, T{ 3.99e-18 } * T{ 1.0e-4 } }, { T{ 192.0e-9 }, T{ 3.60e-18 } * T{ 1.0e-4 } },
      { T{ 194.0e-9 }, T{ 3.51e-18 } * T{ 1.0e-4 } }, { T{ 196.0e-9 }, T{ 3.66e-18 } * T{ 1.0e-4 } },
      { T{ 198.0e-9 }, T{ 3.93e-18 } * T{ 1.0e-4 } }, { T{ 200.0e-9 }, T{ 4.16e-18 } * T{ 1.0e-4 } },
      { T{ 202.0e-9 }, T{ 4.33e-18 } * T{ 1.0e-4 } }, { T{ 204.0e-9 }, T{ 4.40e-18 } * T{ 1.0e-4 } },
      { T{ 206.0e-9 }, T{ 4.45e-18 } * T{ 1.0e-4 } }, { T{ 208.0e-9 }, T{ 4.51e-18 } * T{ 1.0e-4 } },
      { T{ 210.0e-9 }, T{ 4.68e-18 } * T{ 1.0e-4 } }
    };
    return [tabulated](const AtmosphericState<ArrayPolicy> &state, Array3D<T> &weights)
    {
      const auto n_wl = weights.Size1();
      const auto n_z = weights.Size2();
      const auto n_col = weights.Size3();
      for (std::size_t wl = 0; wl < n_wl; ++wl)
      {
        const T lambda_m = state.wavelength_grid_.mid_points_(wl, 0);
        const T wl_nm = lambda_m * T{ 1.0e9 };
        const bool in_band = (wl_nm > T{ 290.0 } && wl_nm < T{ 340.0 });
        for (std::size_t z = 0; z < n_z; ++z)
        {
          for (std::size_t col = 0; col < n_col; ++col)
          {
            const T temperature = state.temperature_(z, col);
            if (in_band && temperature > T{ 210.0 } && temperature < T{ 300.0 })
            {
              const T term1 = (T{ 0.06183 } - (T{ 0.000241 } * wl_nm)) * (T{ 273.0 } - temperature);
              const T term2 = T{ 2.376 } + (T{ 0.14757 } * wl_nm);
              weights(wl, z, col) = std::exp(term1 - term2) * T{ 1.0e-4 };  // cm^2 -> m^2
            }
            else
            {
              weights(wl, z, col) = detail::tabulated_lookup<T>(tabulated, lambda_m);
            }
          }
        }
      }
    };
  }

  /// @brief HNO3 -> OH + NO2 cross-section (m^2): sigma0(l) * exp(sigma1(l)*(T-298)).
  ///
  /// Base sigma0 and exponent coefficient sigma1 are tabulated per wavelength
  /// (from cross_section.hno3.nc). Composed as a base cross-section multiplied
  /// by an exponential temperature scaling.
  template<typename ArrayPolicy = Array3D<double>>
  auto hno3() -> TransformFunc<ArrayPolicy>
  {
    using T = typename ArrayPolicy::value_type;
    auto sigma0 = detail::scaled_array<T>(
        { 1.36e-17, 1.225e-17, 1.095e-17, 9.4e-18, 7.7e-18, 5.88e-18, 4.47e-18, 3.28e-18, 2.31e-18, 1.56e-18,
          1.04e-18 },
        1.0e-4);
    auto sigma1 = detail::first_order_exponent<T>(
        { 0.0017, 0.0017, 0.0017, 0.0017, 0.00165, 0.00166, 0.00169, 0.00174, 0.00177, 0.00185, 0.00197 });
    return tuvx::multiply<ArrayPolicy>(
        tuvx::from_data<ArrayPolicy>(sigma0), tuvx::exponential_scaling<ArrayPolicy>(sigma1, T{ 298.0 }));
  }

  /// @brief RONO2 (alkyl nitrate) cross-section (m^2): sigma0(l) * exp(sigma1(l)*(T-298)).
  template<typename ArrayPolicy = Array3D<double>>
  auto rono2() -> TransformFunc<ArrayPolicy>
  {
    using T = typename ArrayPolicy::value_type;
    auto sigma0 = detail::scaled_array<T>({ 4.65e-21, 3.82e-21, 3.01e-21, 2.42e-21 }, 1.0e-4);
    auto sigma1 = detail::first_order_exponent<T>({ 0.00536, 0.00561, 0.00589, 0.00623 });
    return tuvx::multiply<ArrayPolicy>(
        tuvx::from_data<ArrayPolicy>(sigma0), tuvx::exponential_scaling<ArrayPolicy>(sigma1, T{ 298.0 }));
  }

  /// @brief CH3ONO2 -> CH3O + NO2 cross-section (m^2): sigma0(l) * exp(sigma1(l)*(T-298)).
  template<typename ArrayPolicy = Array3D<double>>
  auto ch3ono2() -> TransformFunc<ArrayPolicy>
  {
    using T = typename ArrayPolicy::value_type;
    auto sigma0 = detail::scaled_array<T>({ 2.98e-21, 2.41e-21, 1.94e-21, 1.54e-21 }, 1.0e-4);
    auto sigma1 = detail::first_order_exponent<T>({ 0.00514, 0.00543, 0.00567, 0.006 });
    return tuvx::multiply<ArrayPolicy>(
        tuvx::from_data<ArrayPolicy>(sigma0), tuvx::exponential_scaling<ArrayPolicy>(sigma1, T{ 298.0 }));
  }

  /// @brief CH2O -> products cross-section (m^2): sigma0(l) + sigma1(l)*(T-298).
  ///
  /// Base value and linear temperature slope are tabulated per wavelength
  /// (from cross_section.ch2o.nc); a straight linear temperature correction.
  template<typename ArrayPolicy = Array3D<double>>
  auto ch2o() -> TransformFunc<ArrayPolicy>
  {
    using T = typename ArrayPolicy::value_type;
    auto base = detail::scaled_array<T>(
        { 4.74e-20, 1.74e-20, 5.56e-20, 1.19e-20, 1.54e-20, 3.86e-20, 9.66e-22, 3.15e-20, 4.36e-21, 3.6e-22,
          6.96e-21, 8.7e-23, 8.8e-23, 6.36e-22 },
        1.0e-4);
    auto slope = detail::scaled_array<T>(
        { 4.173e-24, 2.36e-24, 2.587e-24, 2.4e-25, 2.213e-24, -2.173e-24, -5.0e-27, 5.52e-24, 2.765e-24, -7.28e-25,
          2.456e-24, 0.0, 0.0, 0.0 },
        1.0e-4);
    return tuvx::linear_correction<ArrayPolicy>(base, slope, T{ 298.0 });
  }

  /// @brief CFC-11 -> products cross-section (m^2): sigma0(l) * exp((l_nm - 184.9) * 1e-4 * (T-298)).
  ///
  /// Only the base cross-section is tabulated (from cross_section.cfc-11.nc);
  /// the temperature dependence is a closed-form function of wavelength.
  template<typename ArrayPolicy = Array3D<double>>
  auto cfc11() -> TransformFunc<ArrayPolicy>
  {
    using T = typename ArrayPolicy::value_type;
    auto sigma0 = detail::scaled_array<T>(
        { 1.78e-18, 1.49e-18, 1.23e-18, 9.9e-19, 8.01e-19, 6.47e-19, 5.08e-19, 3.88e-19, 2.93e-19, 2.12e-19,
          1.54e-19 },
        1.0e-4);
    return tuvx::multiply<ArrayPolicy>(
        tuvx::from_data<ArrayPolicy>(sigma0),
        tuvx::bounded_analytic<ArrayPolicy>(
            [](T lambda_m, T temperature) -> T
            {
              const T wl = lambda_m * T{ 1.0e9 };  // nm
              return std::exp(((wl - T{ 184.9 }) * T{ 1.0e-4 }) * (temperature - T{ 298.0 }));
            }));
  }

}  // namespace tuvx::fixed_configuration
