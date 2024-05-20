// Copyright (C) 2023-2024 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
//
/// @file Grid dimensions
#pragma once

#include <string>

#include <tuvx/util/array2d.hpp>

namespace tuvx {

  /// Grid dimensions
  ///
  /// When grid dimensions vary by column, the grid is represented by
  /// two 2D arrays: one for the mid-points and one for the edges,
  /// both indexed by [vertical level][column]
  ///
  /// When grid dimensions are constant, the size of the column
  /// dimension is 1.
  template <typename ArrayPolicy = Array2D<double>>
  class Grid {
  public:
    /// Grid centers
    ArrayPolicy mid_points_;
    /// Grid edges
    ArrayPolicy edges_;

    Grid() = delete;

    /// Constructor for a grid with dimensions that vary by column.
    Grid(const std::string& units, size_t number_of_columns, size_t number_of_grid_cells)
      : units_(units),
        mid_points_(number_of_grid_cells, number_of_columns),
        edges_(number_of_grid_cells+1, number_of_columns),
        is_constant_(false) {}

    /// Constructor for a grid with dimensions that are constant.
    Grid(const std::string& units, size_t number_of_grid_cells)
      : units_(units),
        mid_points_(number_of_grid_cells, 1),
        edges_(number_of_grid_cells+1, 1),
        is_constant_(true) {}

    /// Number of columns
    size_t number_of_columns() const {
      return mid_points_.Size2();
    }

    /// Number of grid cells
    size_t number_of_cells() const {
      return mid_points_.Size1();
    }

    /// Number of grid edges
    size_t number_of_edges() const {
      return edges_.Size1();
    }

    /// Units of the grid.
    std::string units() const {
      return units_;
    }

    /// Check if the grid dimensions are constant across columns.
    bool is_constant() const {
      return is_constant_;
    }
  private:
    /// Units of the grid.
    std::string units_;
    /// True if the grid dimensions are constant across columns.
    bool is_constant_;
  };

} // namespace tuvx