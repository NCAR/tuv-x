#pragma once

#include <iostream>
#include <random>
#include <vector>

namespace tuvx {
namespace linalg {

// tridiag matrix
template <typename T> struct trid_mat {
  int size;
  std::vector<T> udiag; // upper diagonal
  std::vector<T> ldiag; // lower diagonal
  std::vector<T> mdiag; // main diagonal
};

template <typename T> void fill_rand_vec(std::vector<T> &x, int size);

template <typename T> void fill_rand_mat(trid_mat<T> &trid_mat, int n);

template <typename T> void print_vec(std::vector<T>);

template <typename T> void print_trid_mat(trid_mat<T> x);

template <typename T>
std::vector<T> tridiag_solve(trid_mat<T> A, std::vector<T> b);

template <typename T> std::vector<T> dot(trid_mat<T>, std::vector<T>);
} // namespace linalg
} // namespace tuvx
#include "linalg.inl"
