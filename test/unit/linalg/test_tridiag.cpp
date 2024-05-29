#include "../../../include/tuvx/linalg/linalg.h"
#include "../../../src/linalg/tridiag.cpp"

#include <random>

typedef tuvx::linalg::trid_mat<double> trid_matd; 
typedef std::vector<double> vecd;

template <typename T>
void fill_rand(std::vector<T> &x, int size);

template <typename T>
void fill_mat(tuvx::linalg::trid_mat<T> &trid_mat, int n); 

template <typename T>
std::vector<T> dot(tuvx::linalg::trid_mat<T>, std::vector<T>);

template <typename T>
void print_vec(std::vector<T>);

int main() {

    trid_matd M;
    vecd b; 
    vecd x; 
    vecd x_sol; 

    fill_mat(M, 5);
    fill_rand(b, 5);
    fill_rand(x, 5);

    x_sol = tuvx::linalg::tridiag_solve<double>(M, b);

    return 0; 
}

template <typename T>
void fill_rand(std::vector<T> &x, int size){    
    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_int_distribution<std::mt19937::result_type> dist6(1,6);

    x = std::vector<T>(size);
    for(int i = 0; i < x.size(); i++) {
        x[i] = (T)dist6(rng);
    }
}


template <typename T>
void fill_mat(tuvx::linalg::trid_mat<T> &trid_mat, int n) {
    fill_rand(trid_mat.mdiag, n);
    fill_rand(trid_mat.ldiag, n);
    fill_rand(trid_mat.udiag, n);
}

template <typename T>
void print_vec(std::vector<T> x){
    for (int i = 0; i < (int)x.size(); i++)
        std::cout<< x.at(i) <<std::endl;
}

template <typename T>
std::vector<T> dot(tuvx::linalg::trid_mat<T> A, std::vector<T> x){
    std::vector<T> v(x.size());
    v[0] = A.mdiag[0] * x[0] + A.udiag[0]*x[1];

    for (int i = 1; i < x.size(); i++) {
        v[i] = A.mdiag[i]*x[i]   + 
               A.udiag[i-1]*x[i+1] + 
               A.ldiag[i-1]*x[i-1];
    }

    v[i] = A.mdiag[i]*x[i] + A.ldiag[i-1]* 
    return v;
}