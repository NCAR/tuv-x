#pragma once

#include <iostream>
#include <random>
#include <vector>

namespace tuvx {
namespace linalg {

/// @brief Tridiagonal Matrix Data Structure.
/// Given a matrix of shape (n, n), the
/// data structure holds three std::vectors
/// storing the diagonal values.
template <typename T> struct TridiagonalMatrix {
  std::size_t size;
  std::vector<T> upper_diagonal_; // upper diagonal
  std::vector<T> lower_diagonal_; // lower diagonal
  std::vector<T> main_diagonal_;  // main diagonal
};

template <typename T> void fill_rand_vec(std::vector<T> &x, std::size_t size);

template <typename T>
void fill_rand_mat(TridiagonalMatrix<T> &A, std::size_t n);

template <typename T> void print_vec(std::vector<T>);

template <typename T> void print_trid_mat(TridiagonalMatrix<T> x);

template <typename T>
std::vector<T> tridiag_solve(TridiagonalMatrix<T> A, std::vector<T> b);

template <typename T> std::vector<T> dot(TridiagonalMatrix<T>, std::vector<T>);

template <typename T> T norm(std::vector<T> x, int order);

} // namespace linalg
} // namespace tuvx
#include "linalg.inl"
