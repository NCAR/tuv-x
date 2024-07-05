#pragma once

#include <algorithm>
#include <functional>
#include <iostream>
#include <random>
#include <vector>

namespace tuvx
{

  /// @brief Tridiagonal Matrix Data Structure.
  template<typename T>
  struct TridiagonalMatrix
  {
    std::size_t size_;
    std::vector<T> upper_diagonal_;  // upper diagonal
    std::vector<T> lower_diagonal_;  // lower diagonal
    std::vector<T> main_diagonal_;   // main diagonal
    /// @brief tridiagonal matrix structure constructor
    TridiagonalMatrix(std::size_t size);
  };

  /// @brief Fills a vector with uniformly distributed random values.
  /// @param x Vector to and fill
  /// @param seed for random number generation
  template<typename T>
  void FillRandom(std::vector<T> &x, const unsigned &seed);

  /// @brief Fills a matrix with uniformly distributed random values.
  /// @param A tridiagonal matrix to allocate and fill
  /// @param seed for random number generation
  /// @param make_diagonally_dominant Flag that is set to true if matrix is diagonally dominant; false otherwise
  ;
  template<typename T>
  void FillRandom(TridiagonalMatrix<T> &A, const unsigned &seed, const bool &make_diagonally_dominant = false);

  /// @brief display contents of std::vector
  /// @param x Vector to print
  template<typename T>
  void Print(const std::vector<T> &x);

  /// @brief displays the data stored inside a tridiagonal matrix.
  /// @param A tridiagonal matrix to print
  template<typename T>
  void Print(const TridiagonalMatrix<T> &A);

  /// @brief Thomas' algorithm for solving tridiagonal linear system (A x = b).
  /// @param A tridiagonal coeffcient matrix
  /// @param b right hand side vector of the tridiagonal system.
  template<typename T>
  void Solve(TridiagonalMatrix<T> &A, std::vector<T> &b);

  /// @brief Specialized dot product function for tridiagonal matrices.
  /// @param A tridiagonal matrix
  /// @param x vector to multiply the matrix with
  template<typename T>
  std::vector<T> Dot(const TridiagonalMatrix<T> &A, const std::vector<T> &b);

  /// @brief Computes the relative error between two vectors.
  /// @param x true solution
  /// @param x_approx approximated solution
  template<typename T>
  T ComputeError(const std::vector<T> &x, const std::vector<T> &x_approx);

}  // namespace tuvx
#include "linear_algebra.inl"
