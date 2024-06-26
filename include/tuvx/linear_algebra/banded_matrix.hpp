#include <vector>

namespace tuvx
{

  template<typename T>
  struct BandedMatrix
  {
    // data fields
    std::size_t size_;
    std::size_t bandwidth_up_;
    std::size_t bandwidth_down_;
    std::vector<std::vector<T>> band;
    // constructor
    BandedMatrix<T>(std::size_t size, std::size_t up, std::size_t down);
  };

  template<typename T>
  std::vector<T> Dot(const BandedMatrix<T>& A, const std::vector<T>&);

}  // namespace tuvx
