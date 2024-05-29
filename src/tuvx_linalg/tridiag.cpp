#include <tuvx/include/linalg/linalg.h>


namespace tuvx{
  namespace linalg{

    template<typename T>
    vec<T>::vec(size_t N) {
      this->size = N;
      this->elements = (T*)malloc(sizeof(T)*N);
    }

    template<typename T>
    vec<T>::~vec() {
      free(this->elements);
    }

    template<typename T>
    trid_mat() {
      this->size(N);
      this->mdiag = (T*)malloc(sizeof(T)*N);
      this->udiag = (T*)malloc(sizeof(T)*N);
      this->ldiag = (T*)malloc(sizeof(T)*N);
    }

    template <typename T>
    vec<T>::~trid_mat() {
      free(this->mdiag);
      free(this->udiag);
      free(this->ldiag);
    }

    /* solve tridiagonal system */
    template <typename T>
    void* tridiag_solve( trid_mat<T> *A,
	                 dvec<T>     *x,
	                 dvec<T>     *b ) { 
	T denom;
	T temp = A.upper[0] / A.main[0];
 	
	// solve the first lowest diag
	x[0]   = b[0] / A->main[0];  
	
	// forward iteration 
	for (int i = 1; i < b->size; i++) {
	  temp  = A->upper[i] / (A->main[i] - A->lower[i]*temp)
	  x[i]  = b[i] - A->lower[i]*x[i-1];  
        }	

	// backward iteration
	for (int i = b->size-1; i = 1; i--) {
	  b[i] = b[i] - temp*b[i+1];
	  temp = (A->main[i-1] - A->lower[i-1]/temp) / A->upper[i-1]; 
	}
    }

  }	
}
