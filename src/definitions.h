#pragma once
#include <Eigen/Core>

//
// Helpers and miscelaneous functions.
//

namespace TC
{

template <typename Scalar>
inline Scalar shearModulus(Scalar E, Scalar nu) {return E / (2 * (1 + nu));}

//Warning: this function is here just for testing purposes. To be removed soon!
template<typename D1, typename D2>
Eigen::Matrix<typename D1::Scalar, D1::ColsAtCompileTime, D2::ColsAtCompileTime>
  pdist2(const Eigen::MatrixBase<D1>& _X, const Eigen::MatrixBase<D2>& _Y){
  Eigen::Ref<const typename D1::PlainObject> X(_X);
  Eigen::Ref<const typename D2::PlainObject> Y(_Y);
  using matrixXt = Eigen::Matrix<typename D1::Scalar, Eigen::Dynamic, Eigen::Dynamic>;
  using retMatrix = Eigen::Matrix<typename D1::Scalar, D1::ColsAtCompileTime, D2::ColsAtCompileTime>;
  retMatrix dists = X.colwise().squaredNorm().transpose() * matrixXt::Ones(1, Y.cols()) +
    matrixXt::Ones(X.cols(), 1) * Y.colwise().squaredNorm() -
    2 * X.transpose() * Y;
  return dists;
  }

//Warning: this function is here just for testing purposes. To be removed soon!
template<typename D1, typename D2>
Eigen::Matrix<typename D1::Scalar, 3, 3>
p3(const Eigen::MatrixBase<D1>& X, const Eigen::MatrixBase<D2>& Y) {
  //Eigen::Ref<const typename D1::PlainObject> X(_X);
  //Eigen::Ref<const typename D2::PlainObject> Y(_Y);
  using Vector3t = Eigen::Matrix<typename D1::Scalar, 3, 1>;
  using retMatrix = Eigen::Matrix<typename D1::Scalar, 3, 3>;
  EIGEN_STATIC_ASSERT_SAME_MATRIX_SIZE(D1, retMatrix);
  EIGEN_STATIC_ASSERT_SAME_MATRIX_SIZE(D2, retMatrix);
  assert(X.derived().cols() == 3 && X.derived().rows() == 3 &&
         Y.derived().cols() == 3 && Y.derived().rows() == 3);
  retMatrix dists = X.colwise().squaredNorm().transpose() * Vector3t::Ones().transpose() +
    Vector3t::Ones() * Y.colwise().squaredNorm() -
    2 * X.transpose() * Y;
  return dists;
  }

}//namespace TC
