// Copyright (C) 2023-2024 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
//
/// @file Profile of a property on a grid.
#pragma once

#include <string>

#include <tuvx/util/array2d.hpp>

namespace tuvx {

  /// Profile of a property on a grid.
  template <typename T=double>
  class Profile {
  public:
    /// Units of the profile.
    std::string units_;
    /// Values at grid mid-points.
    Array2D<T> mid_point_values_;
    /// Values at grid edges.
    Array2D<T> edge_values_;

    /// Constructor
    Profile(size_t number_of_columns, const Grid& grid)
      : mid_point_values_(number_of_columns, grid.mid_points_.Size1()),
        edge_values_(number_of_columns, grid.edges_.Size1()) {}
  };

} // namespace tuvx