#pragma once
#include <Eigen/Core>
#include "definitions.h"
#include "materialBase.h"

namespace TC
{

template <typename _Scalar, int _Dim>
class AnisotropicLinear {
public:
  typedef typename _Scalar Scalar;
  enum {
    Dim = _Dim,
    StiffDim = _Dim == 3 ? 6 : 3
    };
  using MatrixType = Eigen::Matrix<Scalar, Dim, Dim>;
  using StiffType = Eigen::Matrix<Scalar, StiffDim, StiffDim>;
protected:
  StiffType m_matrix;
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  /** Constructs and initializes the material from a Dim^2 matrix */
  template<typename OtherDerived>
  AnisotropicLinear(const Eigen::MatrixBase<OtherDerived>& other)
    {
    EIGEN_STATIC_ASSERT((Eigen::internal::is_same<Scalar, typename OtherDerived::Scalar>::value),
      YOU_MIXED_DIFFERENT_NUMERIC_TYPES__YOU_NEED_TO_USE_THE_CAST_METHOD_OF_MATRIXBASE_TO_CAST_NUMERIC_TYPES_EXPLICITLY);
    EIGEN_STATIC_ASSERT_SAME_MATRIX_SIZE(StiffType, OtherDerived);
    m_matrix = other.derived();
    }
    
  template<typename Derived>
  MatrixType Stress(const Eigen::MatrixBase<Derived>& strain) const;
  StiffType Stiffness() const {return m_matrix;}
  };

/*template<typename _Scalar, int _Dim>
struct traits<IsotropicLinear<_Scalar, _Dim> >
{
typedef _Scalar Scalar;
enum {
Dim = _Dim,
StiffDim = _Dim = 3 ? 6 : 3
};
};*/

template <typename _Scalar, int _Dim>
template <typename Derived>
typename AnisotropicLinear<_Scalar, _Dim>::MatrixType
AnisotropicLinear<_Scalar, _Dim>::Stress(const Eigen::MatrixBase<Derived>& strain) const
  {
  //I know this is bad. I'll improve it later.
  EIGEN_STATIC_ASSERT_SAME_MATRIX_SIZE(Derived, MatrixType);
  Eigen::Matrix<Scalar, StiffDim,1> strainVec;
  strainVec.head<Dim>() = strain.diagonal();
  strainVec.segment<Dim-1>(Dim) = strain.diagonal<1>();
  if(Dim == 3) strainVec(5) = strain(0,2);
  MatrixType S;
  Eigen::Matrix<Scalar, StiffDim, 1> stressVec = Stiffness() * strainVec; 
  S.diagonal() = stressVec.head<Dim>();
  S.diagonal<1>() = stressVec.segment<Dim - 1>(Dim);
  if (Dim == 3) S(0, 2) = stressVec(5);
  S.triangularView<Eigen::StrictlyLower>() = S.triangularView<Eigen::StrictlyUpper>().transpose();
  return S;// ;
  }

} //namespace TC
#pragma once
