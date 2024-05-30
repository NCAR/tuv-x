#include "../../../include/tuvx/linalg/linalg.h"
#include "../../../src/linalg/linalg.cpp"
#include <random>

typedef tuvx::linalg::trid_mat<double> trid_matd;
typedef std::vector<double> vecd;

template <typename T> void fill_rand(std::vector<T> &x, int size);
template <typename T> void fill_mat(tuvx::linalg::trid_mat<T> &trid_mat, int n);
template <typename T> void print_vec(std::vector<T>);
template <typename T> void print_trid_mat(tuvx::linalg::trid_mat<T> x);

int main() {
  trid_matd M;
  vecd b;
  vecd x = {1, 6, 2, 3, 5};
  M.udiag = {3, 1, 5, 3};
  M.mdiag = {3, 2, 6, 2, 4};
  M.ldiag = {1, 3, 6, 5};
  M.size = 5;
  vecd x_sol;

  print_trid_mat(M);
  print_vec(x);
  b = dot(M, x);
  print_vec(tuvx::linalg::tridiag_solve(M, b));

  /*
  fill_mat(M, 5);
  fill_rand(b, 5);
  fill_rand(x, 5);

  print_trid_mat(M);
  print_vec(x);
  x = tuvx::linalg::dot(M, x);
  x_sol = tuvx::linalg::tridiag_solve<double>(M, b);

  print_vec(x);
  */

  return 0;
}

template <typename T> void fill_rand(std::vector<T> &x, int size) {
  std::random_device dev;
  std::mt19937 rng(dev());
  std::uniform_int_distribution<std::mt19937::result_type> dist6(1, 6);
  x = std::vector<T>(size);
  for (int i = 0; i < x.size(); i++) {
    x[i] = (T)dist6(rng);
  }
}

template <typename T> void fill_mat(tuvx::linalg::trid_mat<T> &tridmat, int n) {
  tridmat.size = n;
  fill_rand(tridmat.mdiag, n);
  fill_rand(tridmat.ldiag, n - 1);
  fill_rand(tridmat.udiag, n - 1);
}

template <typename T> void print_vec(std::vector<T> x) {
  std::cout << std::endl;
  for (int i = 0; i < (int)x.size(); i++) {
    std::cout << x.at(i) << std::endl;
  }
  std::cout << std::endl;
}

template <typename T> void print_trid_mat(tuvx::linalg::trid_mat<T> x) {
  print_vec(x.udiag);
  print_vec(x.mdiag);
  print_vec(x.ldiag);
  std::cout << "----" << std::endl;
}
