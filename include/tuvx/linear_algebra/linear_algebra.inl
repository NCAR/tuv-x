
namespace tuvx
{

  template<typename T>
  inline TridiagonalMatrix<T>::TridiagonalMatrix(std::size_t size)
  {
    this->size_ = size;
    this->main_diagonal_ = std::vector<T>(size);
    this->upper_diagonal_ = std::vector<T>(size - 1);
    this->lower_diagonal_ = std::vector<T>(size - 1);
  }

  template<typename T>
  inline std::vector<T> Dot(const TridiagonalMatrix<T> &A, const std::vector<T> &x)
  {
    std::size_t size = x.size();
    std::vector<T> v(size);
    v[0] = A.main_diagonal_[0] * x[0] + A.upper_diagonal_[0] * x[1];

    std::size_t i = 0;
    for (i = 1; i < size - 1; i++)
    {
      v[i] = A.main_diagonal_[i] * x[i] + A.upper_diagonal_[i] * x[i + 1] + A.lower_diagonal_[i - 1] * x[i - 1];
    }
    v[i] = A.main_diagonal_[i] * x[i] + A.lower_diagonal_[i - 1] * x[i - 1];
    return v;
  }

  template<typename T>
  inline void Solve(TridiagonalMatrix<T> &A, std::vector<T> &b)
  {
    T temp;
    std::size_t N = b.size();
    // forward pass
    for (std::size_t i = 1; i < N; i++)
    {
      temp = A.lower_diagonal_[i - 1] / A.main_diagonal_[i - 1];
      A.main_diagonal_[i] -= temp * A.upper_diagonal_[i - 1];
      b[i] -= temp * b[i - 1];
    }
    // back substitution
    b[N - 1] = b[N - 1] / A.main_diagonal_[N - 1];
    for (std::size_t i = N - 2;; i--)
    {
      b[i] = (b[i] - A.upper_diagonal_[i] * b[i + 1]) / A.main_diagonal_[i];
      if (i == 0)
      {
        break;
      }
    }
  }

  template<typename T>
  inline void FillRandom(std::vector<T> &x)
  {
    std::random_device dev;
    std::mt19937 rng(dev());
    std::uniform_int_distribution<std::mt19937::result_type> dist6(1, 6);
    for (std::size_t i = 0; i < x.size(); i++)
    {
      x[i] = (T)dist6(rng) + 1;
    }
  }

  template<typename T>
  inline void FillRandom(TridiagonalMatrix<T> &A, bool diagonally_dominant)
  {
    FillRandom<T>(A.main_diagonal_);
    FillRandom<T>(A.lower_diagonal_);
    FillRandom<T>(A.upper_diagonal_);

    if (diagonally_dominant)
    {
      // make diagonally dominant (diagonal value greater than sum of its row)
      std::size_t i = 0;
      A.main_diagonal_[i] += A.upper_diagonal_[i];
      for (i = 1; i < A.size_ - 1; i++)
      {
        A.main_diagonal_[i] += A.lower_diagonal_[i - 1] + A.upper_diagonal_[i];
      }
      A.main_diagonal_[i] += A.lower_diagonal_[i - 1];
    }
  }

  template<typename T>
  inline void Print(const std::vector<T> &x)
  {
    std::cout << std::endl;
    for (std::size_t i = 0; i < x.size(); i++)
    {
      std::cout << x.at(i) << std::endl;
    }
    std::cout << std::endl;
  }

  template<typename T>
  T ComputeError(const std::vector<T> &x, const std::vector<T> &x_approx)
  {
    T error = 0;
    for (std::size_t i = 0; i < x.size(); i++)
    {
      error += (abs(x[i] - x_approx[i]) / std::max(x[i], x_approx[i])) / (T)x.size();
    }
    return error;
  }

  template<typename T>
  inline void Print(const TridiagonalMatrix<T> &x)
  {
    std::cout << "----" << std::endl;
    Print(x.upper_diagonal_);
    Print(x.main_diagonal_);
    Print(x.lower_diagonal_);
    std::cout << "----" << std::endl;
  }

}  // namespace tuvx
