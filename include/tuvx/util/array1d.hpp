// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <cstddef>
#include <vector>

namespace tuvx
{

  /// @brief 1D array backed by contiguous storage.
  template<typename T = double>
  class Array1D
  {
   public:
    Array1D() = default;

    /// @brief Construct an array of the given size.
    /// @param size Number of elements.
    explicit Array1D(std::size_t size)
        : data_(size)
    {
    }

    /// @brief Construct an array of the given size filled with a value.
    /// @param size Number of elements.
    /// @param fill_value Value to fill the array with.
    Array1D(std::size_t size, T fill_value)
        : data_(size, fill_value)
    {
    }

    T &operator[](std::size_t i)
    {
      return data_[i];
    }

    const T &operator[](std::size_t i) const
    {
      return data_[i];
    }

    /// @brief Number of elements.
    [[nodiscard]] std::size_t Size() const
    {
      return data_.size();
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

    [[nodiscard]] std::vector<T> &AsVector()
    {
      return data_;
    }
    [[nodiscard]] const std::vector<T> &AsVector() const
    {
      return data_;
    }

   private:
    std::vector<T> data_;
  };

}  // namespace tuvx
