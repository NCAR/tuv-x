
namespace tuvx {

template <typename T>
inline std::vector<T> Dot(TridiagonalMatrix<T> A, std::vector<T> x) {
  std::size_t size = x.size();
  std::vector<T> v(size);
  v[0] = A.main_diagonal_[0] * x[0] + A.upper_diagonal_[0] * x[1];

  for (std::size_t i = 1; i < size - 1; i++) {
    v[i] = A.main_diagonal_[i] * x[i] + A.upper_diagonal_[i] * x[i + 1] +
           A.lower_diagonal_[i - 1] * x[i - 1];
  }
  v[i] = A.main_diagonal_[i] * x[i] + A.lower_diagonal_[i - 1] * x[i - 1];
  return v;
}

template <typename T>
inline std::vector<T> TridiagonalSolve(TridiagonalMatrix<T> A,
                                       std::vector<T> b) {
  T temp;
  std::size_t N = b.size();
  std::vector<T> x(N);
  // forward pass
  for (std::size_t i = 1; i < N; i++) {
    temp = A.lower_diagonal_[i - 1] / A.main_diagonal_[i - 1];
    A.main_diagonal_[i] -= temp * A.upper_diagonal_[i - 1];
    b[i] -= temp * b[i - 1];
  }
  // back substitution
  x[N - 1] = b[N - 1] / A.main_diagonal_[N - 1];
  for (int i = N - 2; i >= 0; i--) {
    x[i] = (b[i] - A.upper_diagonal_[i] * x[i + 1]) / A.main_diagonal_[i];
  }
  return x;
}

template <typename T>
inline void FillRandomVector(std::vector<T> &x, std::size_t size) {
  std::random_device dev;
  std::mt19937 rng(dev());
  std::uniform_int_distribution<std::mt19937::result_type> dist6(1, 6);
  x = std::vector<T>(size);
  for (int i = 0; i < x.size(); i++) {
    x[i] = (T)dist6(rng) + 1;
  }
}

template <typename T>
inline void FillRandomMatrix(TridiagonalMatrix<T> &tridmat, std::size_t n) {
  tridmat.size_ = n;
  FillRandomVector<T>(tridmat.main_diagonal_, n);
  FillRandomVector<T>(tridmat.lower_diagonal_, n - 1);
  FillRandomVector<T>(tridmat.upper_diagonal_, n - 1);
  for (auto i = 0; i < n; i++) {
    tridmat.main_diagonal_[i] = std::pow(tridmat.main_diagonal_[i], 3);
  }
}

template <typename T> inline void Norm(std::vector<T> x, int order) {
  T nrm = (T)0;
  for (auto i = 0; i < x.size(); i++) {
    nrm += std::pow(x[i], order);
  }
  return (1.0 / nrm) * std::pow(nrm, 1.0 / order);
}

template <typename T> inline void PrintVector(std::vector<T> x) {
  std::cout << std::endl;
  for (int i = 0; i < (int)x.size(); i++) {
    std::cout << x.at(i) << std::endl;
  }
  std::cout << std::endl;
}

template <typename T>
inline void PrintTridiagonalMatrix(TridiagonalMatrix<T> x) {
  print_vec(x.upper_diagonal_);
  print_vec(x.main_diagonal_);
  print_vec(x.lower_diagonal_);
  std::cout << "----" << std::endl;
}

} // namespace tuvx
