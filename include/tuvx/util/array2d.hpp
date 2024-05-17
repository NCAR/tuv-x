// Copyright (C) 2023-2024 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
//
/// @file 2D array with row-major storage.
#pragma once

#include <iterator>
#include <vector>

namespace tuvx {

  /// 2D array with row-major storage.
  template <typename T=double>
  class Array2D {
    public:
      Array2D(size_t dim1, size_t dim2)
        : dim1_(dim1), dim2_(dim2), data_(dim1 * dim2) {}

      T &operator()(size_t i, size_t j) { return data_[index(i, j)]; }

      const T &operator()(size_t i, size_t j) const {
        return data_[index(i, j)];
      }

      size_t Size1() const { return dim1_; }
      size_t Size2() const { return dim2_; }

      std::iterator begin() { return data_.begin(); }
      std::iterator end() { return data_.end(); }

    private:
      size_t index(size_t i, size_t j) const { return i * dim2_ + j; }

      size_t dim1_, dim2_;
      std::vector<T> data_;
  };
} // namespace tuvx