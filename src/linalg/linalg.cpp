#include "../../include/tuvx/linalg/linalg.h"

namespace tuvx {
namespace linalg {

template <typename T>
std::vector<T> dot(tuvx::linalg::trid_mat<T> A, std::vector<T> x) {
  int size = x.size();
  std::vector<T> v(size);
  v[0] = A.mdiag[0] * x[0] + A.udiag[0] * x[1];

  int i = 0;
  for (i = 1; i < x.size() - 1; i++) {
    v[i] =
        A.mdiag[i] * x[i] + A.udiag[i] * x[i + 1] + A.ldiag[i - 1] * x[i - 1];
  }
  v[i] = A.mdiag[i] * x[i] + A.ldiag[i - 1] * x[i - 1];
  return v;
}

/* solve tridiagonal system */
template <typename T>
std::vector<T> tridiag_solve(trid_mat<T> A, std::vector<T> b) {
  T temp = A.udiag[0] / A.mdiag[0];
  std::vector<T> x(b.size());

  // solve the first lowest diag
  x[0] = b[0] / A.mdiag[0];

  // forward iteration
  for (int i = 1; i < b.size(); i++) {
    temp = A.udiag[i - 1] / (A.mdiag[i] - A.ldiag[i - 1] * temp);
    x[i] =
        (b[i] - A.ldiag[i] * x[i - 1]) / (A.mdiag[i] - A.ldiag[i - 1] * temp);
  }

  // backward iteration
  for (int i = b.size(); i > 0; i--) {
    b[i] = b[i] - temp * b[i + 1];
    temp = (A.mdiag[i] - (A.ldiag[i - 1] / temp)) / A.ldiag[i - 1];
  }
  return x;
}
/*
    cp(1) = c(1) / b(1)
    u(1) = r(1) / b(1)
    do i = 2, size( b )
      denom = 1.0 / ( b(i) - a(i) * cp(i-1) )
      cp(i) = c(i) * denom
      u(i) = ( r(i) - a(i) * u(i-1) ) * denom
    end do
    do i = size( b ) - 1, 1, -1
      u(i) = u(i) - cp(i) * u(i+1)
    end do

 */
} // namespace linalg
} // namespace tuvx
