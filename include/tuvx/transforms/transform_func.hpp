// Copyright (C) 2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
//
// TransformFunc - the open, user-facing API for calculating weight arrays.
#pragma once

#include <tuvx/radiative_transfer/atmospheric_state.hpp>
#include <tuvx/util/array3d.hpp>

#include <functional>

namespace tuvx
{

  /// @brief Callable that calculates a transform weight array.
  ///
  /// A TransformFunc fills `weights[wavelength x height x column]` given the
  /// current atmospheric state.  Cross-sections, quantum yields, and spectral
  /// weights are all TransformFuncs.  Any callable matching this signature is a
  /// valid transform - raw lambdas, function objects, or free functions - no
  /// factory wrapper is required.
  ///
  /// The weight array is *written* by the transform (it does not accumulate into
  /// an existing value).  Combinators that need to combine two transforms
  /// allocate a temporary internally.
  ///
  /// Example - user-written lambda:
  /// @code
  ///   tuvx::TransformFunc<> rayleigh = [](const auto& state, auto& weights) {
  ///       weights.ForEachRow(
  ///           [](double& w, const double& wl) {
  ///               constexpr double A = 1.0455996e-30;  // m2 nm4
  ///               w = A / (wl * wl * wl * wl);
  ///           },
  ///           weights.GetColumnView(0),
  ///           state.wavelength_grid_.mid_points_);
  ///   };
  /// @endcode
  ///
  /// @tparam ArrayPolicy  3D array policy (default: Array3D<double>).
  template<typename ArrayPolicy = Array3D<double>>
  using TransformFunc = std::function<void(
      const AtmosphericState<ArrayPolicy> &state,
      Array3D<typename ArrayPolicy::value_type> &weights)>;

}  // namespace tuvx
