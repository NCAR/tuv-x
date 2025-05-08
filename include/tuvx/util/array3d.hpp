// Copyright (C) 2023-2025 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <cstddef>
#include <iterator>
#include <vector>

namespace tuvx
{

  /// @brief 3D array with row-major storage.
  template<typename T = double>
  class Array3D
  {
   public:
    Array3D() = default;

    Array3D(std::size_t dim1, std::size_t dim2, std::size_t dim3)
        : dim1_(dim1),
          dim2_(dim2),
          dim3_(dim3),
          data_(dim1 * dim2 * dim3)
    {
    }

    T &operator()(std::size_t i, std::size_t j, std::size_t k)
    {
      return data_[index(i, j, k)];
    }

    const T &operator()(std::size_t i, std::size_t j, std::size_t k) const
    {
      return data_[index(i, j, k)];
    }

    std::size_t Size1() const
    {
      return dim1_;
    }

    std::size_t Size2() const
    {
      return dim2_;
    }

    std::size_t Size3() const
    {
      return dim3_;
    }

    typename std::vector<T>::iterator begin()
    {
      return data_.begin();
    }

    typename std::vector<T>::iterator end()
    {
      return data_.end();
    }

    typename std::vector<T>::const_iterator begin() const
    {
      return data_.begin();
    }

    typename std::vector<T>::const_iterator end() const
    {
      return data_.end();
    }

    std::vector<T> &AsVector()
    {
      return data_;
    }
    const std::vector<T> &AsVector() const
    {
      return data_;
    }

   private:
    std::size_t index(std::size_t i, std::size_t j, std::size_t k) const
    {
      return i * dim2_ * dim3_ + j * dim3_ + k;
    }

    std::size_t dim1_, dim2_, dim3_;
    std::vector<T> data_;
  };

}  // namespace tuvx