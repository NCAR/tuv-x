// Copyright (C) 2023-2024 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
//
/// @file Profile of a property on a grid.
#pragma once

#include <string>

#include <tuvx/util/array2d.hpp>
#include <tuvx/grid.hpp>

namespace tuvx {

  /// Profile of a property on a grid.
  template <typename ArrayPolicy = Array2D<double>>
  class Profile {
  public:
    /// Values at grid mid-points.
    ArrayPolicy mid_point_values_;
    /// Values at grid edges.
    ArrayPolicy edge_values_;

    Profile() = delete;

    /// Constructor
    template<typename U>
    Profile(const std::string& units, size_t number_of_columns, const Grid<U>& grid)
      : units_(units),
        mid_point_values_(grid.mid_points_.Size1(), number_of_columns),
        edge_values_(grid.edges_.Size1(), number_of_columns) {}

    /// Units of the profile.
    std::string units() const {
      return units_;
    }

  private:
    /// Units of the profile.
    std::string units_;

  };

} // namespace tuvx