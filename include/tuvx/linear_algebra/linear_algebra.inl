
#include <random>
namespace tuvx
{

  template<typename T>
  inline TridiagonalMatrix<T>::TridiagonalMatrix(std::size_t size)
      : size_(size),
        main_diagonal_(size),
        upper_diagonal_(size - 1),
        lower_diagonal_(size - 1)
  {
  }

  template<typename T>
  inline Array1D<T> Dot(const TridiagonalMatrix<T> &A, const Array1D<T> &x)
  {
    std::size_t size = x.Size();
    Array1D<T> v(size);
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
  inline void Solve(TridiagonalMatrix<T> &A, Array1D<T> &b)
  {
    T temp;
    std::size_t n = b.Size();
    // forward pass
    for (std::size_t i = 1; i < n; i++)
    {
      temp = A.lower_diagonal_[i - 1] / A.main_diagonal_[i - 1];
      A.main_diagonal_[i] -= temp * A.upper_diagonal_[i - 1];
      b[i] -= temp * b[i - 1];
    }
    // back substitution
    b[n - 1] = b[n - 1] / A.main_diagonal_[n - 1];
    for (std::size_t i = n - 2;; i--)
    {
      b[i] = (b[i] - A.upper_diagonal_[i] * b[i + 1]) / A.main_diagonal_[i];
      if (i == 0)
      {
        break;
      }
    }
  }

  template<typename T>
  inline void FillRandom(Array1D<T> &x, const unsigned &seed)
  {
    std::mt19937 random_device(seed);
    std::normal_distribution<double> distribution(5.0, 1.0);
    for (std::size_t i = 0; i < x.Size(); i++)
    {
      x[i] = (T)distribution(random_device);
    }
  }

  template<typename T>
  inline void FillRandom(TridiagonalMatrix<T> &A, const unsigned &seed, const bool &make_diagonally_dominant)
  {
    FillRandom<T>(A.main_diagonal_, seed);
    FillRandom<T>(A.lower_diagonal_, seed);
    FillRandom<T>(A.upper_diagonal_, seed);

    if (make_diagonally_dominant)
    {
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
  inline void Print(const Array1D<T> &x)
  {
    std::cout << std::endl;
    for (std::size_t i = 0; i < x.Size(); i++)
    {
      std::cout << x[i] << std::endl;
    }
    std::cout << std::endl;
  }

  template<typename T>
  T ComputeError(const Array1D<T> &x, const Array1D<T> &x_approx)
  {
    T error = 0;
    for (std::size_t i = 0; i < x.Size(); i++)
    {
      error += (std::abs(x[i] - x_approx[i]) / std::max(x[i], x_approx[i])) / (T)x.Size();
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
