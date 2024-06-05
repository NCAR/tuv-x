
namespace tuvx {
namespace linalg {

template <typename T>
inline std::vector<T> dot(trid_mat<T> A, std::vector<T> x) {
  std::size_t size = x.size();
  std::vector<T> v(size);
  v[0] = A.mdiag[0] * x[0] + A.udiag[0] * x[1];

  std::size_t i = 0;
  for (i = 1; i < size - 1; i++) {
    v[i] =
        A.mdiag[i] * x[i] + A.udiag[i] * x[i + 1] + A.ldiag[i - 1] * x[i - 1];
  }
  v[i] = A.mdiag[i] * x[i] + A.ldiag[i - 1] * x[i - 1];
  return v;
}

template <typename T>
inline std::vector<T> tridiag_solve(trid_mat<T> A, std::vector<T> b) {
  T temp;
  std::size_t N = b.size();
  std::vector<T> x(N);
  // forward pass
  for (int i = 1; i < N; i++) {
    temp = A.ldiag[i - 1] / A.mdiag[i - 1];
    A.mdiag[i] -= temp * A.udiag[i - 1];
    b[i] -= temp * b[i - 1];
  }
  // back substitution
  x[N - 1] = b[N - 1] / A.mdiag[N - 1];
  for (int i = N - 2; i >= 0; i--) {
    x[i] = (b[i] - A.udiag[i] * x[i + 1]) / A.mdiag[i];
  }
  return x;
}

template <typename T>
inline void fill_rand_vec(std::vector<T> &x, std::size_t size) {
  std::random_device dev;
  std::mt19937 rng(dev());
  std::uniform_int_distribution<std::mt19937::result_type> dist6(1, 6);
  x = std::vector<T>(size);
  for (int i = 0; i < x.size(); i++) {
    x[i] = (T)dist6(rng) + 1;
  }
}

template <typename T>
inline void fill_rand_mat(trid_mat<T> &tridmat, std::size_t n) {
  tridmat.size = n;
  fill_rand_vec(tridmat.mdiag, n);
  fill_rand_vec(tridmat.ldiag, n - 1);
  fill_rand_vec(tridmat.udiag, n - 1);
  for (int i = 0; i < n; i++) {
    tridmat.mdiag[i] = std::pow(tridmat.mdiag[i], 3);
  }
}

template <typename T> inline void norm(std::vector<T> x, int order) {
  T nrm = (T)0;
  for (int i = 0; i < x.size(); i++) {
    nrm += std::pow(x[i], order);
  }
  return (1.0 / nrm) * std::pow(nrm, 1.0 / order);
}

template <typename T> inline void print_vec(std::vector<T> x) {
  std::cout << std::endl;
  for (int i = 0; i < (int)x.size(); i++) {
    std::cout << x.at(i) << std::endl;
  }
  std::cout << std::endl;
}

template <typename T> inline void print_trid_mat(trid_mat<T> x) {
  print_vec(x.udiag);
  print_vec(x.mdiag);
  print_vec(x.ldiag);
  std::cout << "----" << std::endl;
}

} // namespace linalg
} // namespace tuvx
