#include "../../include/tuvx/linalg/linalg.h"

namespace tuvx
{
  namespace linalg
  {

    /* solve tridiagonal system */
    template <typename T>
    std::vector<T> tridiag_solve(trid_mat<T> A, std::vector<T> b) {
      T denom;
      T temp = A.udiag[0] / A.mdiag[0];
      std::vector<T> x(b.size()); 

      // solve the first lowest diag
      x[0] = b[0] / A.mdiag[0];

      // forward iteration
      for (int i = 1; i < b.size(); i++) {
        temp = A.udiag[i] / (A.mdiag[i] - A.ldiag[i] * temp);
        x[i] = b[i] - A.ldiag[i] * x[i - 1];
      }

      // backward iteration
      for (int i = b.size() - 1; i>0; i--) {
        b[i] = b[i] - temp * b[i + 1];
        temp = (A.mdiag[i - 1] - A.ldiag[i - 1] / temp) / A.ldiag[i - 1];
      }
      return x;
    }
  }
}
