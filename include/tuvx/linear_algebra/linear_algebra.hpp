#pragma once

#include <tuvx/util/array1d.hpp>

#include <algorithm>
#include <functional>
#include <iostream>
#include <random>

namespace tuvx
{

  /// @brief Tridiagonal Matrix Data Structure.
  template<typename T>
  struct TridiagonalMatrix
  {
    std::size_t size_;
    Array1D<T> upper_diagonal_;  // upper diagonal
    Array1D<T> lower_diagonal_;  // lower diagonal
    Array1D<T> main_diagonal_;   // main diagonal
    /// @brief Initializes the three internal vectors based on
    /// size. Main diagonal will be of size $n$, while the upper/lower
    /// diagonals will be of size $n-1$
    /// @param size Size of the matrix to be initialized
    TridiagonalMatrix(std::size_t size);
  };

  /// @brief Fills an Array1D with random values
  /// @param x Array to fill
  /// @param seed Seed for random number generation
  template<typename T>
  void FillRandom(Array1D<T> &x, const unsigned &seed);

  /// @brief Fills a matrix with uniformly distributed random values.
  /// @param A Tridiagonal matrix to fill
  /// @param seed Seed for random number generation
  /// @param make_diagonally_dominant Make the tridiagonal matrix diagonally dominant
  template<typename T>
  void FillRandom(TridiagonalMatrix<T> &A, const unsigned &seed, const bool &make_diagonally_dominant = false);

  /// @brief Displays the data stored inside an Array1D
  /// @param x Array to print
  template<typename T>
  void Print(const Array1D<T> &x);

  /// @brief Displays the data stored inside a tridiagonal matrix.
  /// @param A Tridiagonal matrix to print
  template<typename T>
  void Print(const TridiagonalMatrix<T> &x);

  /// @brief Thomas' algorithm for solving tridiagonal linear system (A x = b).
  /// Solution is stored in b.
  /// @param A Tridiagonal coefficient matrix
  /// @param b Right hand side vector of the tridiagonal system.
  template<typename T>
  void Solve(TridiagonalMatrix<T> &A, Array1D<T> &b);

  /// @brief Specialized dot product function for tridiagonal matrices.
  /// @param A Tridiagonal matrix
  /// @param x Array to multiply the matrix with
  /// @returns Dot product between A and x
  template<typename T>
  Array1D<T> Dot(const TridiagonalMatrix<T> &A, const Array1D<T> &x);

  /// @brief Computes the relative error between two arrays.
  /// @param x True solution
  /// @param x_approx Approximated solution
  template<typename T>
  T ComputeError(const Array1D<T> &x, const Array1D<T> &x_approx);

}  // namespace tuvx
#include "linear_algebra.inl"
