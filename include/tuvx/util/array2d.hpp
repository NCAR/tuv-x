// Copyright (C) 2023-2024 National Science Foundation-National Center for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <cstddef>
#include <vector>

namespace tuvx
{

  /// @brief 2D array with row-major storage.
  template<typename T = double>
  class Array2D
  {
   public:
    Array2D() = default;

    Array2D(std::size_t dim1, std::size_t dim2)
        : dim1_(dim1),
          dim2_(dim2),
          data_(dim1 * dim2)
    {
    }

    T &operator()(std::size_t i, std::size_t j)
    {
      return data_[index(i, j)];
    }

    const T &operator()(std::size_t i, std::size_t j) const
    {
      return data_[index(i, j)];
    }

    std::size_t Size1() const
    {
      return dim1_;
    }
    std::size_t Size2() const
    {
      return dim2_;
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
    std::size_t index(std::size_t i, std::size_t j) const
    {
      return i * dim2_ + j;
    }

    std::size_t dim1_, dim2_;
    std::vector<T> data_;
  };
}  // namespace tuvx