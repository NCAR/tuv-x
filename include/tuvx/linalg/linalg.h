#pragma once

#include <stdlib.h>

namespace tuvx{
namespace linalg{

  // dense vector
  template <typename T> 
  class vec{
    private:
      size_t size; 
      T* elements; 
    public:
      vec(size_t size);
      ~vec();
  };

  // tridiag matrix
  template <typename T>
  class trid_mat {
    private:
      size_t size;
      T* udiag; // upper diagonal 
      T* ldiag; // lower diagonal
      T* mdiag; // main diagonal
  }

  template <typename T>
  void tridiag_solve(trid_mat<T> A, vec<T> x, vec<T> b);

} 
}
