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
  Scalar _E;
  Scalar _nu;
public:
  IsotropicLinear(Scalar E, Scalar nu): _E(E), _nu(nu)
    {};
  template<typename Derived>
  MatrixType Stress(const Eigen::MatrixBase<Derived>& strain) const;
  StiffType Stiffness() const;
  Scalar E() const { return _E; }
  Scalar nu() const { return _nu; }
  Scalar G() const {return .5*E() / (1 + nu());}
  };

template <typename _Scalar, int _Dim>
typename IsotropicLinear<_Scalar, _Dim>::StiffType 
    IsotropicLinear<_Scalar, _Dim>::Stiffness() const
  {
  StiffType K = StiffType::Zero();
  Scalar ee = E()/((1+nu())*(1-2*nu()));
  K.template topLeftCorner<Dim, Dim>() = MatrixType::Constant(nu()*ee) + (1-2*nu())*ee * MatrixType::Identity();
  K.diagonal().template tail<StiffDim - Dim>() =  Eigen::Matrix<Scalar, StiffDim - Dim, 1>::Constant(G());
  return K;
  }

template <typename _Scalar, int _Dim>
template <typename Derived>
typename IsotropicLinear<_Scalar, _Dim>::MatrixType 
    IsotropicLinear<_Scalar, _Dim>::Stress(const Eigen::MatrixBase<Derived>& strain) const
  {
  EIGEN_STATIC_ASSERT_SAME_MATRIX_SIZE(Derived, MatrixType);
  MatrixType S;
  Scalar ee = E()/((1.+nu())*(1.-2.*nu()));
  Eigen::Matrix<_Scalar, Dim, 1> diag = (MatrixType::Constant(nu()*ee) + (1-2*nu())*ee * MatrixType::Identity()) * strain.derived().diagonal();
  S = (MatrixType::Constant(G()) - G() * MatrixType::Identity()).cwiseProduct(strain.derived());
  S += diag.asDiagonal();
  return S;// ;
  }

} //namespace TC
