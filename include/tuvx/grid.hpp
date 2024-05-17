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
  template <typename T=double>
  class Grid {
  public:
    /// Units of the grid.
    std::string units_;
    /// Grid centers
    Array2D<T> mid_points_;
    /// Grid edges
    Array2D<T> edges_;

    /// Constructor for a grid with dimensions that vary by column.
    Grid(size_t number_of_columns, size_t number_of_grid_cells)
      : mid_points_(number_of_grid_cells, number_of_columns),
        edges_(number_of_grid_cells+1, number_of_columns),
        is_constant_(false) {}

    /// Constructor for a grid with dimensions that are constant.
    Grid(size_t number_of_grid_cells)
      : mid_points_(number_of_grid_cells, 1),
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

    /// Check if the grid dimensions are constant across columns.
    bool is_constant() const {
      return is_constant_;
    }
  private:
    bool is_constant_;
  };

} // namespace tuvx