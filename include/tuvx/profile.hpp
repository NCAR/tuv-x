// Copyright (C) 2023-2025 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <tuvx/grid.hpp>
#include <tuvx/util/array2d.hpp>

#include <cstddef>
#include <string>

namespace tuvx
{

  /// @brief Profile of a property on a grid.
  template<typename ArrayPolicy = Array2D<double>>
  class Profile
  {
   public:
    /// Values at grid mid-points.
    ArrayPolicy mid_point_values_;
    /// Values at grid edges.
    ArrayPolicy edge_values_;

    Profile() = default;

    /// @brief Constructor of a profile
    /// @param units Units of the profile.
    /// @param number_of_columns Number of columns.
    /// @param grid Grid that the profile is defined on.
    template<typename GridPolicy>
    Profile(const std::string& units, std::size_t number_of_columns, const GridPolicy& grid)
        : units_(units),
          mid_point_values_(grid.mid_points_.Size1(), number_of_columns),
          edge_values_(grid.edges_.Size1(), number_of_columns)
    {
    }

    /// Units of the profile.
    std::string Units() const
    {
      return units_;
    }

   private:
    /// Units of the profile.
    std::string units_;
  };

}  // namespace tuvx