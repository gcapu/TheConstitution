#pragma once
#include <Eigen/Core>
#include "definitions.h"
#include "materialBase.h"

namespace TC
{
template <typename _Scalar, int _Dim>
class IsotropicLinear {
public:
  typedef _Scalar Scalar;

  enum { 
      Dim = _Dim,
      StiffDim = _Dim==3?6:3
  };
  using MatrixType = Eigen::Matrix<Scalar, Dim, Dim>;
  using StiffType = Eigen::Matrix<Scalar, StiffDim, StiffDim>;
protected:
  Scalar _lambda;
  Scalar _mu;
  Scalar _density;
public:
  IsotropicLinear(Scalar E, Scalar nu, Scalar density = 1.): 
      _lambda(nu*E/(1.+nu)/(1.-2.*nu)), _mu(E/2./(1+nu)), _density(density) {};
  template<typename Derived>
  inline MatrixType Stress(const Eigen::MatrixBase<Derived>& strain) const;
  inline StiffType Stiffness() const;
  Scalar E() const { return _mu*(3*_lambda + 2*_mu)/(_lambda + _mu); }
  Scalar nu() const { return _lambda/2./(_lambda+_mu); }
  Scalar mu() const {return _mu;}
  Scalar lambda() const {return _lambda;}
  Scalar density() const {return _density;}
  };

template <typename _Scalar, int _Dim>
typename IsotropicLinear<_Scalar, _Dim>::StiffType 
    IsotropicLinear<_Scalar, _Dim>::Stiffness() const
  {
  StiffType K = mu()* StiffType::Identity();
  K.template topLeftCorner<Dim, Dim>() = MatrixType::Constant(lambda()) + 2.*mu() * MatrixType::Identity();
  return K;
  }

template <typename _Scalar, int _Dim>
template <typename Derived>
typename IsotropicLinear<_Scalar, _Dim>::MatrixType 
    IsotropicLinear<_Scalar, _Dim>::Stress(const Eigen::MatrixBase<Derived>& strain) const
  {
  EIGEN_STATIC_ASSERT_SAME_MATRIX_SIZE(Derived, MatrixType);
  return lambda() * MatrixType::Identity() * strain.trace() + 2 * mu() * strain;
  }

} //namespace TC
