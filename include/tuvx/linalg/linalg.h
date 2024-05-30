#pragma once

#include <iostream>
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

template <typename T>
std::vector<T> tridiag_solve(trid_mat<T> A, std::vector<T> b);

template <typename T>
std::vector<T> dot(tuvx::linalg::trid_mat<T>, std::vector<T>);
} // namespace linalg
} // namespace tuvx
