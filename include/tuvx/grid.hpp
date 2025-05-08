// Copyright (C) 2023-2025 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <tuvx/util/array2d.hpp>

#include <string>

namespace tuvx
{

  /// @brief Grid dimensions
  ///
  /// When grid dimensions vary by column, the grid is represented by
  /// two 2D arrays: one for the mid-points and one for the edges,
  /// both indexed by [vertical level][column]
  ///
  /// When grid dimensions are constant, the size of the column
  /// dimension is 1.
  template<typename ArrayPolicy = Array2D<double>>
  class Grid
  {
   public:
    /// Grid centers
    ArrayPolicy mid_points_;
    /// Grid edges
    ArrayPolicy edges_;

    Grid() = default;

    /// @brief Constructor for a grid with dimensions that vary by column.
    /// @param units Units of the grid.
    /// @param number_of_columns Number of columns.
    /// @param number_of_sections Number of grid sections per column.
    Grid(const std::string& units, std::size_t number_of_columns, std::size_t number_of_sections)
        : units_(units),
          mid_points_(number_of_sections, number_of_columns),
          edges_(number_of_sections + 1, number_of_columns),
          is_constant_(false)
    {
    }

    /// @brief Constructor for a grid with dimensions that are constant.
    /// @param units Units of the grid.
    /// @param number_of_sections Number of grid sections per column.
    Grid(const std::string& units, std::size_t number_of_sections)
        : units_(units),
          mid_points_(number_of_sections, 1),
          edges_(number_of_sections + 1, 1),
          is_constant_(true)
    {
    }

    /// @brief Number of columns
    ///
    /// When the grid dimensions are constant, the number of columns is 1.
    std::size_t NumberOfColumns() const
    {
      return mid_points_.Size2();
    }

    /// @brief Number of grid sections
    ///
    /// The number of grid sections refers to the discreet intervals along
    /// the grid for a single column. The total number of values stored by
    /// the grid would then be the number of grid sections times the number
    /// of columns.
    std::size_t NumberOfSections() const
    {
      return mid_points_.Size1();
    }

    /// @brief Number of grid edges
    ///
    /// Grid edges are the boundaries of the grid sections. The number of
    /// grid edges is one more than the number of grid sections.
    std::size_t NumberOfEdges() const
    {
      return edges_.Size1();
    }

    /// @brief Units of the grid.
    std::string Units() const
    {
      return units_;
    }

    /// @brief Check if the grid dimensions are constant across columns.
    bool IsConstant() const
    {
      return is_constant_;
    }

   private:
    /// Units of the grid.
    std::string units_;
    /// True if the grid dimensions are constant across columns.
    bool is_constant_;
  };

}  // namespace tuvx