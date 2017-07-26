#include <Eigen/Core>

template <typename Derived>
typename Eigen::EigenBase<Derived>::Index cols(const Eigen::EigenBase<Derived>& a)
  {
  return a.cols();
  }

template <typename Derived1, typename Derived2>
inline Eigen::Matrix& matrixstuff(const Eigen::MatrixBase<Derived1>& a, const Eigen::MatrixBase<Derived2>& b)
  {
  return a*b;
  }
