// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
//
// Convenience factory functions that return TransformFunc callables.
//
// Each factory captures its parameters and returns a lambda matching the
// TransformFunc signature.  Users are free to write raw TransformFunc lambdas
// directly - these factories are helpers, not required entry points.
#pragma once

#include <tuvx/transforms/transform_func.hpp>
#include <tuvx/util/array1d.hpp>
#include <tuvx/util/array2d.hpp>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <concepts>
#include <cstddef>

namespace tuvx
{

  // ---------------------------------------------------------------------------
  // constant

  /// @brief Returns a TransformFunc that sets every weight to @p value.
  ///
  /// Useful for flat quantum yields or unit weights.
  ///
  /// @param value  Weight value broadcast over all [wavelength x height x column].
  /// @return TransformFunc that fills @p weights with @p value.
  template<typename ArrayPolicy = Array3D<double>>
  auto constant(typename ArrayPolicy::value_type value) -> TransformFunc<ArrayPolicy>
  {
    return [value](const AtmosphericState<ArrayPolicy> & /*state*/, Array3D<typename ArrayPolicy::value_type> &weights)
    {
      for (auto &w : weights)
      {
        w = value;
      }
    };
  }

  // ---------------------------------------------------------------------------
  // wrap_analytic
  //
  // Two overloads selected by concept:
  //   wrap_analytic(f)  where f(lambda) -> T     -- wavelength-only formula
  //   wrap_analytic(f)  where f(lambda, T) -> T  -- wavelength + temperature formula

  /// @brief Concept for a wavelength-only analytic formula.
  template<typename F, typename T>
  concept WavelengthFormula = requires(F f, T wl) {
    { f(wl) } -> std::convertible_to<T>;
  };

  /// @brief Concept for a wavelength-and-temperature analytic formula.
  template<typename F, typename T>
  concept WavelengthTemperatureFormula = requires(F f, T wl, T temp) {
    { f(wl, temp) } -> std::convertible_to<T>;
  };

  /// @brief Wraps a simple analytic formula f(lambda) -> T into a full TransformFunc.
  ///
  /// The result is broadcast over all height levels and columns.
  ///
  /// Example:
  /// @code
  ///   auto xs = tuvx::wrap_analytic([](double wl) {
  ///       return 1.0e-22 * std::exp(-(wl - 300e-9) * (wl - 300e-9) / (10e-9 * 10e-9));
  ///   });
  /// @endcode
  ///
  /// @tparam ArrayPolicy  Storage policy.
  /// @param  f  Callable returning a weight given wavelength (m).
  template<typename ArrayPolicy = Array3D<double>, typename F>
    requires WavelengthFormula<F, typename ArrayPolicy::value_type> &&
                 (!WavelengthTemperatureFormula<F, typename ArrayPolicy::value_type>)
  auto wrap_analytic(F f) -> TransformFunc<ArrayPolicy>
  {
    return [f = std::move(f)](const AtmosphericState<ArrayPolicy> &state, Array3D<typename ArrayPolicy::value_type> &weights)
    {
      const auto n_wl = weights.Size1();
      const auto n_z = weights.Size2();
      const auto n_col = weights.Size3();
      for (std::size_t wl = 0; wl < n_wl; ++wl)
      {
        const auto lambda = state.wavelength_grid_.mid_points_(wl, 0);
        const auto w = f(lambda);
        for (std::size_t z = 0; z < n_z; ++z)
        {
          for (std::size_t col = 0; col < n_col; ++col)
          {
            weights(wl, z, col) = w;
          }
        }
      }
    };
  }

  /// @brief Wraps an analytic formula f(lambda, T) -> T into a full TransformFunc.
  ///
  /// Temperature is taken from AtmosphericState::temperature_ at each (height, column).
  ///
  /// @tparam ArrayPolicy  Storage policy.
  /// @param  f  Callable returning a weight given wavelength (m) and temperature (K).
  template<typename ArrayPolicy = Array3D<double>, typename F>
    requires WavelengthTemperatureFormula<F, typename ArrayPolicy::value_type>
  auto wrap_analytic(F f) -> TransformFunc<ArrayPolicy>
  {
    return [f = std::move(f)](const AtmosphericState<ArrayPolicy> &state, Array3D<typename ArrayPolicy::value_type> &weights)
    {
      const auto n_wl = weights.Size1();
      const auto n_z = weights.Size2();
      const auto n_col = weights.Size3();
      for (std::size_t wl = 0; wl < n_wl; ++wl)
      {
        const auto lambda = state.wavelength_grid_.mid_points_(wl, 0);
        for (std::size_t z = 0; z < n_z; ++z)
        {
          for (std::size_t col = 0; col < n_col; ++col)
          {
            weights(wl, z, col) = f(lambda, state.temperature_(z, col));
          }
        }
      }
    };
  }

  // ---------------------------------------------------------------------------
  // from_data

  /// @brief Returns a TransformFunc from a wavelength-indexed table of values.
  ///
  /// @p model_values must already be interpolated onto the model wavelength grid
  /// (length == number of wavelength bins the transform will be called with).
  /// Values are broadcast uniformly over all height levels and columns.
  ///
  /// @param model_values  Pre-interpolated values indexed by wavelength bin [n_wavelengths].
  template<typename ArrayPolicy = Array3D<double>>
  auto from_data(Array1D<typename ArrayPolicy::value_type> model_values) -> TransformFunc<ArrayPolicy>
  {
    return [vals = std::move(model_values)](
               const AtmosphericState<ArrayPolicy> & /*state*/, Array3D<typename ArrayPolicy::value_type> &weights)
    {
      assert(vals.Size() == weights.Size1());
      const auto n_z = weights.Size2();
      const auto n_col = weights.Size3();
      for (std::size_t wl = 0; wl < vals.Size(); ++wl)
      {
        for (std::size_t z = 0; z < n_z; ++z)
        {
          for (std::size_t col = 0; col < n_col; ++col)
          {
            weights(wl, z, col) = vals[wl];
          }
        }
      }
    };
  }

  // ---------------------------------------------------------------------------
  // temperature_interpolation

  /// @brief Returns a TransformFunc using T-dependent tabulated data.
  ///
  /// The @p interp callable encapsulates a reference cross-section table and its
  /// wavelength/temperature grids.  At each call it receives the model wavelength
  /// bin index and the local temperature and returns the interpolated value.
  ///
  /// Example - linear interpolation in T:
  /// @code
  ///   auto xs = tuvx::temperature_interpolation(
  ///       [ref_temps, xs_at_wl](std::size_t wl_idx, double T) {
  ///           return interpolate_linear(ref_temps, xs_at_wl[wl_idx], T);
  ///       });
  /// @endcode
  ///
  /// @tparam ArrayPolicy  Storage policy.
  /// @param  interp  Callable: (wl_idx, temperature_K) -> value.
  template<typename ArrayPolicy = Array3D<double>, typename Interp>
  auto temperature_interpolation(Interp interp) -> TransformFunc<ArrayPolicy>
  {
    return [interp = std::move(interp)](
               const AtmosphericState<ArrayPolicy> &state, Array3D<typename ArrayPolicy::value_type> &weights)
    {
      const auto n_wl = weights.Size1();
      const auto n_z = weights.Size2();
      const auto n_col = weights.Size3();
      for (std::size_t wl = 0; wl < n_wl; ++wl)
      {
        for (std::size_t z = 0; z < n_z; ++z)
        {
          for (std::size_t col = 0; col < n_col; ++col)
          {
            weights(wl, z, col) = interp(wl, state.temperature_(z, col));
          }
        }
      }
    };
  }

  // ---------------------------------------------------------------------------
  // polynomial_scaling
  //
  // Computes sum_n coeffs(wl, n) * (T - T_ref)^n

  /// @brief Returns a TransformFunc for polynomial temperature scaling.
  ///
  /// At each (wavelength, height, column), computes:
  /// \f[ w = \sum_n c_n(\lambda) \left( T - T_\mathrm{ref} \right)^n \f]
  ///
  /// The zero-th order coefficient coeffs(wl, 0) is the base value
  /// \f$ \sigma_0(\lambda) \f$.
  ///
  /// @param coeffs  Polynomial coefficients [n_wavelengths x (max_order+1)].
  ///                coeffs(wl, n) is the n-th order coefficient for wavelength wl.
  /// @param t_ref   Reference temperature (K).
  template<typename ArrayPolicy = Array3D<double>>
  auto polynomial_scaling(Array2D<typename ArrayPolicy::value_type> coeffs, typename ArrayPolicy::value_type t_ref)
      -> TransformFunc<ArrayPolicy>
  {
    return [coeffs = std::move(coeffs), t_ref](
               const AtmosphericState<ArrayPolicy> &state, Array3D<typename ArrayPolicy::value_type> &weights)
    {
      assert(coeffs.Size1() == weights.Size1());
      const auto n_wl = weights.Size1();
      const auto n_z = weights.Size2();
      const auto n_col = weights.Size3();
      const auto poly_order = coeffs.Size2();
      for (std::size_t wl = 0; wl < n_wl; ++wl)
      {
        for (std::size_t z = 0; z < n_z; ++z)
        {
          for (std::size_t col = 0; col < n_col; ++col)
          {
            const auto delta_t = state.temperature_(z, col) - t_ref;
            typename ArrayPolicy::value_type val = 0;
            typename ArrayPolicy::value_type power = 1;
            for (std::size_t n = 0; n < poly_order; ++n)
            {
              val += coeffs(wl, n) * power;
              power *= delta_t;
            }
            weights(wl, z, col) = val;
          }
        }
      }
    };
  }

  // ---------------------------------------------------------------------------
  // exponential_scaling
  //
  // Computes exp(sum_n coeffs(wl, n) * (T - T_ref)^n)

  /// @brief Returns a TransformFunc for exponential temperature scaling.
  ///
  /// At each (wavelength, height, column), computes:
  /// \f[ w = \exp\!\left( \sum_n c_n(\lambda) \left( T - T_\mathrm{ref} \right)^n \right) \f]
  ///
  /// @param coeffs  Polynomial exponent coefficients [n_wavelengths x (max_order+1)].
  /// @param t_ref   Reference temperature (K).
  template<typename ArrayPolicy = Array3D<double>>
  auto exponential_scaling(Array2D<typename ArrayPolicy::value_type> coeffs, typename ArrayPolicy::value_type t_ref)
      -> TransformFunc<ArrayPolicy>
  {
    return [coeffs = std::move(coeffs), t_ref](
               const AtmosphericState<ArrayPolicy> &state, Array3D<typename ArrayPolicy::value_type> &weights)
    {
      assert(coeffs.Size1() == weights.Size1());
      const auto n_wl = weights.Size1();
      const auto n_z = weights.Size2();
      const auto n_col = weights.Size3();
      const auto poly_order = coeffs.Size2();
      for (std::size_t wl = 0; wl < n_wl; ++wl)
      {
        for (std::size_t z = 0; z < n_z; ++z)
        {
          for (std::size_t col = 0; col < n_col; ++col)
          {
            const auto delta_t = state.temperature_(z, col) - t_ref;
            typename ArrayPolicy::value_type exponent = 0;
            typename ArrayPolicy::value_type power = 1;
            for (std::size_t n = 0; n < poly_order; ++n)
            {
              exponent += coeffs(wl, n) * power;
              power *= delta_t;
            }
            weights(wl, z, col) = std::exp(exponent);
          }
        }
      }
    };
  }

  // ---------------------------------------------------------------------------
  // linear_correction
  //
  // Computes base(wl) + slope(wl) * (T - T_ref)

  /// @brief Returns a TransformFunc for a linear temperature correction.
  ///
  /// At each (wavelength, height, column), computes:
  /// \f[ w = \sigma_0(\lambda) + \beta(\lambda) \left( T - T_\mathrm{ref} \right) \f]
  ///
  /// @param base   Base values at reference temperature [n_wavelengths].
  /// @param slope  Temperature slope per wavelength bin [n_wavelengths] (value K-1).
  /// @param t_ref  Reference temperature (K).
  template<typename ArrayPolicy = Array3D<double>>
  auto linear_correction(
      Array1D<typename ArrayPolicy::value_type> base,
      Array1D<typename ArrayPolicy::value_type> slope,
      typename ArrayPolicy::value_type t_ref) -> TransformFunc<ArrayPolicy>
  {
    return [base = std::move(base), slope = std::move(slope), t_ref](
               const AtmosphericState<ArrayPolicy> &state, Array3D<typename ArrayPolicy::value_type> &weights)
    {
      assert(base.Size() == weights.Size1());
      assert(slope.Size() == weights.Size1());
      const auto n_wl = weights.Size1();
      const auto n_z = weights.Size2();
      const auto n_col = weights.Size3();
      for (std::size_t wl = 0; wl < n_wl; ++wl)
      {
        for (std::size_t z = 0; z < n_z; ++z)
        {
          for (std::size_t col = 0; col < n_col; ++col)
          {
            weights(wl, z, col) = base[wl] + (slope[wl] * (state.temperature_(z, col) - t_ref));
          }
        }
      }
    };
  }

  // ---------------------------------------------------------------------------
  // stern_volmer
  //
  // Computes phi = phi0 / (1 + k * M * phi0)
  // where M is air number density (mol m-3) from AtmosphericState::air_density_.

  /// @brief Returns a TransformFunc for Stern-Volmer quenching quantum yield.
  ///
  /// At each (wavelength, height, column), computes:
  /// \f[ w = \frac{\phi_0}{1 + k M \phi_0} \f]
  ///
  /// @param phi0  Base (pressure-independent) quantum yield.
  /// @param k     Quenching rate coefficient (m3 mol-1).
  template<typename ArrayPolicy = Array3D<double>>
  auto stern_volmer(typename ArrayPolicy::value_type phi0, typename ArrayPolicy::value_type k) -> TransformFunc<ArrayPolicy>
  {
    return [phi0, k](const AtmosphericState<ArrayPolicy> &state, Array3D<typename ArrayPolicy::value_type> &weights)
    {
      const auto n_wl = weights.Size1();
      const auto n_z = weights.Size2();
      const auto n_col = weights.Size3();
      for (std::size_t wl = 0; wl < n_wl; ++wl)
      {
        for (std::size_t z = 0; z < n_z; ++z)
        {
          for (std::size_t col = 0; col < n_col; ++col)
          {
            const auto air = state.air_density_(z, col);
            weights(wl, z, col) = phi0 / (typename ArrayPolicy::value_type{ 1 } + (k * air * phi0));
          }
        }
      }
    };
  }

  // ---------------------------------------------------------------------------
  // parameterized
  //
  // Generic escape hatch: the user supplies a callable that, given
  // (wl_idx, z_idx, col_idx, temperature_K, air_density_mol_m-3),
  // returns the weight value.

  /// @brief Returns a TransformFunc from a generic per-element calculator.
  ///
  /// The @p calc callable receives (wl_idx, z_idx, col_idx, temperature_K,
  /// air_density_mol_m3) and returns the weight at that position.
  ///
  /// @tparam ArrayPolicy  Storage policy.
  /// @param  calc  Per-element weight calculator.
  template<typename ArrayPolicy = Array3D<double>, typename Calc>
  auto parameterized(Calc calc) -> TransformFunc<ArrayPolicy>
  {
    return [calc = std::move(calc)](
               const AtmosphericState<ArrayPolicy> &state, Array3D<typename ArrayPolicy::value_type> &weights)
    {
      const auto n_wl = weights.Size1();
      const auto n_z = weights.Size2();
      const auto n_col = weights.Size3();
      for (std::size_t wl = 0; wl < n_wl; ++wl)
      {
        for (std::size_t z = 0; z < n_z; ++z)
        {
          for (std::size_t col = 0; col < n_col; ++col)
          {
            weights(wl, z, col) = calc(wl, z, col, state.temperature_(z, col), state.air_density_(z, col));
          }
        }
      }
    };
  }

}  // namespace tuvx
