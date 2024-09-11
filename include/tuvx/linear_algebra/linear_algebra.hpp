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
    /// @brief Initializes the three internal vectors based on
    /// size. Main diagonal will be of size $n$, while the upper/lower
    /// diagonals will be of size $n-1$
    /// @param size Size of the matrix to be initialized
    TridiagonalMatrix(std::size_t size);
  };

  /// @brief Fills a std::vector with random values
  /// @param x Vector to allocate and fill
  /// @param seed Seed for random number generation
  template<typename T>
  void FillRandom(std::vector<T> &x, const unsigned &seed);

  /// @brief Fills a matrix with uniformly distributed random values.
  /// @param A Tridiagonal matrix to allocate and fill
  /// @param seed Seed for random number generation
  /// @param make_diagonally_dominant Make the tridiagonal matrix diagonally dominant (diagonal value is greater than sum of
  /// the corrosponding row values)
  template<typename T>
  void FillRandom(TridiagonalMatrix<T> &A, const unsigned &seed, const bool &make_diagonally_dominant = false);

  /// @brief Displays the data stored inside a std::vector
  /// @param x Vector to print
  template<typename T>
  void Print(const std::vector<T> &x);

  /// @brief Displays the data stored inside a tridiagonal matrix.
  /// For now, this function simply calls printvec() on the three diagonals.
  /// @param A Tridiagonal matrix to print
  template<typename T>
  void Print(const TridiagonalMatrix<T> &A);

  /// @brief Thomas' algorithm for solving tridiagonal linear system (A x = b).
  /// store solution in b.
  /// @param A Tridiagonal coeffcient matrix
  /// @param b Right hand side vector of the tridiagonal system.
  template<typename T>
  void Solve(TridiagonalMatrix<T> &A, std::vector<T> &b);

  /// @brief Specialized dot product function for tridiagonal matrices.
  /// @param A Tridiagonal matrix
  /// @param x Vector to multiply the matrix with
  /// @returns Dot product between A and x
  template<typename T>
  std::vector<T> Dot(const TridiagonalMatrix<T> &A, const std::vector<T> &b);

  /// @brief Computes the relative error between two vectors. Used for computing
  /// approximation errors.
  /// @param x True solution
  /// @param x_approx Approximated solution
  template<typename T>
  T ComputeError(const std::vector<T> &x, const std::vector<T> &x_approx);

}  // namespace tuvx
#include "linear_algebra.inl"
