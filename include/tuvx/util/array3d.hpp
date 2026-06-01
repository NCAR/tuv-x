// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <cassert>
#include <cstddef>
#include <functional>
#include <vector>

namespace tuvx
{

  /// @brief 3D array with row-major storage [dim1 × dim2 × columns].
  ///
  /// The column dimension is always the last (innermost) dimension.
  /// Bulk operations access data through column views rather than direct
  /// element indexing so that the layout can be changed by policy subclasses
  /// without modifying algorithm code.
  ///
  /// For transforms, the typical layout is [wavelength × height × column].
  template<typename T = double>
  class Array3D
  {
   public:
    using value_type = T;

    // -------------------------------------------------------------------------
    // Column views — lightweight non-owning proxies into one column's data.
    //
    // Views hold a non-owning pointer to the parent array; they are only valid
    // for the lifetime of that array.

    /// @brief Mutable view of one column across all (dim1, dim2) positions.
    class ColumnView
    {
     public:
      ColumnView(Array3D<T> *arr, std::size_t col)
          : arr_(arr),
            col_(col)
      {
      }
      T &operator()(std::size_t i, std::size_t j)
      {
        return (*arr_)(i, j, col_);
      }
      const T &operator()(std::size_t i, std::size_t j) const
      {
        return (*arr_)(i, j, col_);
      }

     private:
      Array3D<T> *arr_ = nullptr;
      std::size_t col_ = 0;
    };

    /// @brief Read-only view of one column across all (dim1, dim2) positions.
    class ConstColumnView
    {
     public:
      ConstColumnView(const Array3D<T> *arr, std::size_t col)
          : arr_(arr),
            col_(col)
      {
      }
      const T &operator()(std::size_t i, std::size_t j) const
      {
        return (*arr_)(i, j, col_);
      }

     private:
      const Array3D<T> *arr_ = nullptr;
      std::size_t col_ = 0;
    };

    /// @brief Temporary per-(dim1,dim2)-pair scalar storage; participates in ForEachRow like a column view.
    class RowVariable
    {
     public:
      RowVariable(std::size_t dim1, std::size_t dim2)
          : data_((dim1 * dim2)),
            dim2_(dim2)
      {
      }
      T &operator()(std::size_t i, std::size_t j)
      {
        return data_[(i * dim2_) + j];
      }
      const T &operator()(std::size_t i, std::size_t j) const
      {
        return data_[(i * dim2_) + j];
      }

     private:
      std::vector<T> data_{};
      std::size_t dim2_ = 0;
    };

    // -------------------------------------------------------------------------
    // Function — reusable, policy-compiled operation over this array type.
    //
    // For the CPU host policy, Function is a lightweight wrapper: it captures
    // the callable and re-applies it on demand. GPU policies pre-compile
    // kernels at construction time using the prototype's dimension metadata.
    //
    // Usage:
    //   auto f = Array3D<double>::Function(
    //       [](auto& arr) {
    //           auto tmp = arr.GetRowVariable();
    //           arr.ForEachRow(..., arr.GetConstColumnView(0), tmp);
    //           arr.ForEachRow(..., tmp, arr.GetColumnView(1));
    //       },
    //       prototype_array);
    //   f(array_a);   // apply to any compatible array

    /// @brief Reusable policy-compiled operation.
    class Function
    {
     public:
      /// @brief Construct a reusable function from a callable and a prototype array.
      /// @param f          Callable that operates on an Array3D via ForEachRow/ColumnView.
      /// @param prototype  Prototype array (reserved for GPU policy pre-compilation).
      template<typename F>
      Function(F f, const Array3D & /*prototype*/)
          : f_(std::move(f))
      {
      }

      /// @brief Apply the compiled function to an array.
      void operator()(Array3D &arr) const
      {
        f_(arr);
      }

     private:
      std::function<void(Array3D<T> &)> f_;
    };

    // -------------------------------------------------------------------------

    Array3D() = default;

    /// @brief Construct an array with the given dimensions.
    /// @param dim1 Size of the first dimension.
    /// @param dim2 Size of the second dimension.
    /// @param dim3 Size of the third (column) dimension.
    Array3D(std::size_t dim1, std::size_t dim2, std::size_t dim3)
        : dim1_(dim1),
          dim2_(dim2),
          dim3_(dim3),
          data_(dim1 * dim2 * dim3)
    {
    }

    T &operator()(std::size_t i, std::size_t j, std::size_t k)
    {
      return data_[(((i * dim2_) + j) * dim3_) + k];
    }

    const T &operator()(std::size_t i, std::size_t j, std::size_t k) const
    {
      return data_[(((i * dim2_) + j) * dim3_) + k];
    }

    /// @brief Size of the first dimension.
    [[nodiscard]] std::size_t Size1() const
    {
      return dim1_;
    }

    /// @brief Size of the second dimension.
    [[nodiscard]] std::size_t Size2() const
    {
      return dim2_;
    }

    /// @brief Size of the third (column) dimension.
    [[nodiscard]] std::size_t Size3() const
    {
      return dim3_;
    }

    // -------------------------------------------------------------------------
    // Bulk operations API

    /// @brief Return a mutable column view for bulk iteration.
    /// @param col Column index [0, Size3()).
    ColumnView GetColumnView(std::size_t col)
    {
      assert(col < dim3_);
      return ColumnView(this, col);
    }

    /// @brief Return a read-only column view for bulk iteration.
    /// @param col Column index [0, Size3()).
    ConstColumnView GetConstColumnView(std::size_t col) const
    {
      assert(col < dim3_);
      return ConstColumnView(this, col);
    }

    /// @brief Return temporary per-(dim1,dim2)-pair storage for intermediate calculations.
    RowVariable GetRowVariable() const
    {
      return RowVariable(dim1_, dim2_);
    }

    /// @brief Iterate over all (dim1, dim2) pairs, invoking @p f with one element per view per pair.
    ///
    /// The lambda receives references to the corresponding element of each
    /// supplied view at each (i, j) index pair.  Views may be ColumnView,
    /// ConstColumnView, RowVariable, Array2D, or any object supporting
    /// operator()(std::size_t, std::size_t).
    ///
    /// Example (scale column 0 into column 1 of a weights array):
    /// @code
    ///   weights.ForEachRow(
    ///       [](const double& a, double& b) { b = a * 2.0; },
    ///       weights.GetConstColumnView(0),
    ///       weights.GetColumnView(1));
    /// @endcode
    ///
    /// Passing an Array2D<T> as a view broadcasts it across both (i, j) dimensions
    /// since Array2D already provides operator()(std::size_t, std::size_t):
    /// @code
    ///   // temperature is Array2D<double>[height × col]; use as view in a 3D ForEachRow
    ///   weights.ForEachRow(
    ///       [](double& w, const double& t) { w *= factor(t); },
    ///       weights.GetColumnView(col),
    ///       temperature);
    /// @endcode
    template<typename Func, typename... Views>
    void ForEachRow(Func &&f, Views &&...views)  // NOLINT(cppcoreguidelines-missing-std-forward)
    {
      for (std::size_t i = 0; i < dim1_; ++i)
      {
        for (std::size_t j = 0; j < dim2_; ++j)
        {
          std::invoke(std::forward<Func>(f), views(i, j)...);
        }
      }
    }

    // -------------------------------------------------------------------------

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

    T *Data()
    {
      return data_.data();
    }
    const T *Data() const
    {
      return data_.data();
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
    std::size_t dim1_ = 0;
    std::size_t dim2_ = 0;
    std::size_t dim3_ = 0;
    std::vector<T> data_{};
  };

}  // namespace tuvx
