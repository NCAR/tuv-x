// Copyright (C) 2023-2024 National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
//
/// @file 3D array with row-major storage.
#pragma once

#include <iterator>
#include <vector>

namespace tuvx {

/// @brief 3D array with row-major storage.
template <typename T = double> class Array3D {
public:
  Array3D() = default;

  Array3D(size_t dim1, size_t dim2, size_t dim3)
      : dim1_(dim1), dim2_(dim2), dim3_(dim3), data_(dim1 * dim2 * dim3) {}

  T &operator()(size_t i, size_t j, size_t k) { return data_[index(i, j, k)]; }

  const T &operator()(size_t i, size_t j, size_t k) const {
    return data_[index(i, j, k)];
  }

  size_t Size1() const { return dim1_; }

  size_t Size2() const { return dim2_; }

  size_t Size3() const { return dim3_; }

  typename std::vector<T>::iterator begin() { return data_.begin(); }

  typename std::vector<T>::iterator end() { return data_.end(); }

  typename std::vector<T>::const_iterator begin() const {
    return data_.begin();
  }

  typename std::vector<T>::const_iterator end() const { return data_.end(); }

  std::vector<T> &AsVector() { return data_; }
  const std::vector<T> &AsVector() const { return data_; }

private:
  size_t index(size_t i, size_t j, size_t k) const {
    return i * dim2_ * dim3_ + j * dim3_ + k;
  }

  size_t dim1_, dim2_, dim3_;
  std::vector<T> data_;
};

} // namespace tuvx