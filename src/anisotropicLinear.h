#pragma once
#include <memory>
#include <Eigen/Core>
#include "utils.h"
namespace TC
{

//Anisotropic linear material stores a matrix and uses it as its stiffness. 
//  If you have multiple copies of the same material it is better to use its
//  copy constructor since it will just copy a pointer to the same matrix.
template <typename _Scalar, int _Dim>
class AnisotropicLinear {
public:
  typedef _Scalar Scalar;
  enum {
    Dim = _Dim,
    StiffDim = _Dim == 3 ? 6 : 3
    };
  using StressType = Eigen::Matrix<Scalar, Dim, Dim>;
  using StiffType = Eigen::Matrix<Scalar, StiffDim, StiffDim>;
protected:
  std::shared_ptr<StiffType> _matrix;
  Scalar _density;
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  // Constructs and initializes the material from a Dim^2 matrix 
  template<typename Derived>
  AnisotropicLinear(const Eigen::MatrixBase<Derived>& stiffness);
  // Copies of the same material point to the same matrix
  AnisotropicLinear(const AnisotropicLinear<Scalar, Dim>& other);
  //Common functions
  Scalar density() const {return _density;}
  template<typename Derived>
  inline StressType Stress(const Eigen::MatrixBase<Derived>& strain) const;
  inline StiffType Stiffness() const {return *_matrix;}
  };

template <typename _Scalar, int _Dim>
template <typename Derived>
AnisotropicLinear<_Scalar, _Dim>::AnisotropicLinear(const Eigen::MatrixBase<Derived>& stiffness)
  {
  EIGEN_STATIC_ASSERT((Eigen::internal::is_same<Scalar, typename Derived::Scalar>::value),
    YOU_MIXED_DIFFERENT_NUMERIC_TYPES__YOU_NEED_TO_USE_THE_CAST_METHOD_OF_MATRIXBASE_TO_CAST_NUMERIC_TYPES_EXPLICITLY);
  EIGEN_STATIC_ASSERT_SAME_MATRIX_SIZE(StiffType, Derived);
  _matrix = std::make_shared<StiffType>(stiffness.derived());
  }
  
template <typename _Scalar, int _Dim>
AnisotropicLinear<_Scalar, _Dim>::AnisotropicLinear(const AnisotropicLinear<Scalar, Dim>& other)
  {
  _matrix = other._matrix;
  _density = other._density;
  }


template <typename _Scalar, int _Dim>
template <typename Derived>
typename AnisotropicLinear<_Scalar, _Dim>::StressType
AnisotropicLinear<_Scalar, _Dim>::Stress(const Eigen::MatrixBase<Derived>& strain) const
  {
  EIGEN_STATIC_ASSERT_SAME_MATRIX_SIZE(Derived, StressType);
  Eigen::Matrix<_Scalar, _Dim, _Dim> strainVal = strain.derived();
  Eigen::Matrix<Scalar, StiffDim,1> strainVec = StrainToVoigt(strainVal);
  Eigen::Matrix<Scalar, StiffDim, 1> stressVec = Stiffness() * strainVec; 
  return VoigtToStress(stressVec);
  }
  
} //namespace TC
