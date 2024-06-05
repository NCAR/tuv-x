#pragma once

#include <iostream>
#include <random>
#include <vector>

namespace tuvx {

/// @struct Tridiagonal Matrix Data Structure.
/// @brief Given a matrix of shape (n, n), the
/// data structure holds three std::vectors
/// storing the diagonal values.
template <typename T> struct TridiagonalMatrix {
  std::size_t size;
  std::vector<T> upper_diagonal_; // upper diagonal
  std::vector<T> lower_diagonal_; // lower diagonal
  std::vector<T> main_diagonal_;  // main diagonal
};

/// @brief Allocates a vector to a given size and fills it with random values
/// @param x Vector to allocate and fill
/// @param size Size of the vector
template <typename T> void FillRandVec(std::vector<T> &x, std::size_t size);

/// @brief Allocates a matrix to a given size and fills it with random values
/// @param A tridiagonal matrix to allocate and fill
/// @param size Size of the vector
template <typename T> void FillRandMat(TridiagonalMatrix<T> &A, std::size_t n);

/// @brief displays the data stored inside a std::vector
/// @param x Vector to print
template <typename T> void PrintVec(std::vector<T>);

/// @brief displays the data stored inside a tridiagonal matrix.
/// For now, this function simply calls printvec() on the three diagonals.
/// @param A tridiagonal matrix to print
template <typename T> void PrintTridMat(TridiagonalMatrix<T> x);

/// @fn Tridiagonal Linear System Solver
/// @brief Thomas' algorithm for solving tridiagonal linear system (A x = b)
/// @param A tridiagonal coeffcient matrix
/// @param b right hand side vector of the tridiagonal system.
/// @returns x solution that solves A x = b.
template <typename T>
std::vector<T> TridiagSolve(TridiagonalMatrix<T> A, std::vector<T> b);

/// @fn Tridiagonal Matrix - vector dot product
/// @brief Specialized dot product function for tridiagonal matrices.
/// @param A tridiagonal matrix
/// @param x vector to multiply the matrix with
/// @returns dot product between A and x
template <typename T> std::vector<T> Dot(TridiagonalMatrix<T>, std::vector<T>);

/// @fn LP norm
/// @brief LP of a vectors. Used for computing errors
/// @param x vector to compute the norm of
/// @param p of the norm
/// @returns $\frac{1}{N}(\sum_{i=1}^{N} |x|^p)^{\frac{1}{p}}$
template <typename T> T Norm(std::vector<T> x, int p);

} // namespace tuvx
#include "linalg.inl"
