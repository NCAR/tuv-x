// Copyright (C) 2023-2026 University Corporation for Atmospheric Research
// SPDX-License-Identifier: Apache-2.0
#pragma once

#include <cassert>
#include <cstddef>
#include <functional>
#include <vector>

namespace tuvx
{

  /// @brief 2D array with row-major storage [rows × columns].
  ///
  /// Dimensions are referred to as "rows" (dim1) and "columns" (dim2).
  /// Bulk operations access data through column views rather than direct
  /// element indexing so that the layout can be changed by policy subclasses
  /// without modifying algorithm code.
  template<typename T = double>
  class Array2D
  {
   public:
    using value_type = T;

    // -------------------------------------------------------------------------
    // Column views — lightweight non-owning proxies into one column's data.
    //
    // Views hold a non-owning pointer to the parent array; they are only valid
    // for the lifetime of that array.

    /// @brief Mutable view of one column across all rows.
    class ColumnView
    {
     public:
      ColumnView(Array2D<T> *arr, std::size_t col)
          : arr_(arr),
            col_(col)
      {
      }
      T &operator[](std::size_t row)
      {
        return (*arr_)(row, col_);
      }
      const T &operator[](std::size_t row) const
      {
        return (*arr_)(row, col_);
      }

     private:
      Array2D<T> *arr_ = nullptr;
      std::size_t col_ = 0;
    };

    /// @brief Read-only view of one column across all rows.
    class ConstColumnView
    {
     public:
      ConstColumnView(const Array2D<T> *arr, std::size_t col)
          : arr_(arr),
            col_(col)
      {
      }
      const T &operator[](std::size_t row) const
      {
        return (*arr_)(row, col_);
      }

     private:
      const Array2D<T> *arr_ = nullptr;
      std::size_t col_ = 0;
    };

    /// @brief Temporary per-row scalar storage; participates in ForEachRow like a column view.
    class RowVariable
    {
     public:
      explicit RowVariable(std::size_t n_rows)
          : data_(n_rows)
      {
      }
      T &operator[](std::size_t row)
      {
        return data_[row];
      }
      const T &operator[](std::size_t row) const
      {
        return data_[row];
      }

     private:
      std::vector<T> data_{};
    };

    // -------------------------------------------------------------------------
    // Function — reusable, policy-compiled operation over this array type.
    //
    // For the CPU host policy, Function is a lightweight wrapper: it captures
    // the callable and re-applies it on demand. GPU policies pre-compile
    // kernels at construction time using the prototype's dimension metadata.
    //
    // Usage:
    //   auto f = Array2D<double>::Function(
    //       [](auto& arr) { arr.ForEachRow(...); },
    //       prototype_array);
    //   f(array_a);   // apply to any compatible array
    //   f(array_b);

    /// @brief Reusable policy-compiled operation.
    class Function
    {
     public:
      /// @brief Construct a reusable function from a callable and a prototype array.
      /// @param f          Callable that operates on an Array2D via ForEachRow/ColumnView.
      /// @param prototype  Prototype array (reserved for GPU policy pre-compilation).
      template<typename F>
      Function(F f, const Array2D & /*prototype*/)
          : f_(std::move(f))
      {
      }

      /// @brief Apply the compiled function to an array.
      void operator()(Array2D &arr) const
      {
        f_(arr);
      }

     private:
      std::function<void(Array2D<T> &)> f_;
    };

    // -------------------------------------------------------------------------

    Array2D() = default;

    /// @brief Construct an array with the given dimensions.
    /// @param dim1 Number of rows.
    /// @param dim2 Number of columns.
    Array2D(std::size_t dim1, std::size_t dim2)
        : dim1_(dim1),
          dim2_(dim2),
          data_(dim1 * dim2)
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

    /// @brief Number of rows.
    [[nodiscard]] std::size_t Size1() const
    {
      return dim1_;
    }

    /// @brief Number of columns.
    [[nodiscard]] std::size_t Size2() const
    {
      return dim2_;
    }

    // -------------------------------------------------------------------------
    // Bulk operations API

    /// @brief Return a mutable column view for bulk iteration.
    /// @param col Column index [0, Size2()).
    ColumnView GetColumnView(std::size_t col)
    {
      assert(col < dim2_);
      return ColumnView(this, col);
    }

    /// @brief Return a read-only column view for bulk iteration.
    /// @param col Column index [0, Size2()).
    ConstColumnView GetConstColumnView(std::size_t col) const
    {
      assert(col < dim2_);
      return ConstColumnView(this, col);
    }

    /// @brief Return temporary per-row storage for intermediate calculations.
    RowVariable GetRowVariable() const
    {
      return RowVariable(dim1_);
    }

    /// @brief Iterate over all rows, invoking @p f with one element per view per row.
    ///
    /// The lambda receives references to the corresponding element of each
    /// supplied view at each row index.  Views may be ColumnView, ConstColumnView,
    /// RowVariable, Array1D, or any object supporting operator[](std::size_t).
    ///
    /// Example (scale column 0 into column 1):
    /// @code
    ///   arr.ForEachRow(
    ///       [](const double& a, double& b) { b = a * 2.0; },
    ///       arr.GetConstColumnView(0),
    ///       arr.GetColumnView(1));
    /// @endcode
    template<typename Func, typename... Views>
    void ForEachRow(Func &&f, Views &&...views)  // NOLINT(cppcoreguidelines-missing-std-forward)
    {
      for (std::size_t i = 0; i < dim1_; ++i)
      {
        std::invoke(std::forward<Func>(f), views[i]...);
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
    std::vector<T> data_{};
  };

}  // namespace tuvx
