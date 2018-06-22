#pragma once
#include <Eigen/Core>
#include "utils.h"
namespace TC
{

//NeoHookean material 
template <typename _Scalar, int _Dim>
class NeoHookean {
public:
  typedef _Scalar Scalar;
  enum {
    Dim = _Dim,
    StiffDim = _Dim == 3 ? 6 : 3
    };
  using StressType = Eigen::Matrix<Scalar, Dim, Dim>;
  using StiffType = Eigen::Matrix<Scalar, StiffDim, StiffDim>;
  struct ResultType {
    StressType S;
    StiffType C;
  };
protected:
  Scalar _lambda;
  Scalar _kappa;
public:
  NeoHookean(Scalar lambda, Scalar kappa): 
      _lambda(lambda), _kappa(kappa) {};
  template<typename Derived>
  inline StressType Stress(const Eigen::MatrixBase<Derived>& F) const;
  template<typename Derived>
  inline ResultType Stiffness(const Eigen::MatrixBase<Derived>& F) const;
  };

template <typename _Scalar, int _Dim>
template <typename Derived>
typename NeoHookean<_Scalar, _Dim>::StressType
NeoHookean<_Scalar, _Dim>::Stress(const Eigen::MatrixBase<Derived>& F) const
  {
  EIGEN_STATIC_ASSERT_SAME_MATRIX_SIZE(Derived, StressType);
  Eigen::Matrix<_Scalar, _Dim, _Dim> strainVal = strain.derived();
  Eigen::Matrix<Scalar, StiffDim,1> strainVec = StrainToVoigt(strainVal);
  Eigen::Matrix<Scalar, StiffDim, 1> stressVec = Stiffness() * strainVec; 
  return VoigtToStress(stressVec);
  }
template <typename _Scalar, int _Dim>
typename NeoHookean<_Scalar, _Dim>::ResultType
NeoHookean<_Scalar, _Dim>::Stiffness(const Eigen::MatrixBase<Derived>& F) const
  {
    static_assert(_Dim==2, "sorry, Neo hookean material is only made for 2d plane strain for now");
    ResultType result;
    //These arrays allow to convert an voigt index to square matrix index
    std::array<int, 3> index1 = {{0,1,0}};
    std::array<int, 3> index2 = {{0,1,1}};
    //C and its inverse
    Eigen::Matrix2d C = F.transpose()*F;
    Eigen::Matrix2d Ci = C.inverse();
    //invariants of C (the second invariant is not necessary)
    double I1 = C.trace()+1; //It's plane strain so the thrid component of C is 1
    double I3 = C.determinant();
    //derivatives of the invariants wrt C in voigt form
    Eigen::Vector3d dI1; dI1 << 1, 1, 0;
    Eigen::Vector3d dI3; dI3 <<I3*Ci(0,0), I3*Ci(1,1), I3*Ci(0,1);
    //second derivatives of the invariants wrt C (the second derivative of I1 is zero)
    Eigen::Matrix3d d2I3; //in voigt notation 
    for(int i = 0; i<3; i++)
      for(int j = 0; j<3; j++)
        {
        int k = index1[i], l = index2 [i];
        int m = index1[j], n = index2 [j];
        d2I3(i,j) = I3*(Ci(k,l)*Ci(m,n)-.5*(Ci(k,m)*Ci(n,l)+Ci(k,n)*Ci(m,l)));
        }
    //reduced invariants of C (the second invariant is also not necessary)
    //double Ib = I1 * pow(I3, -1./3.); //only its derivative is used
    double J = sqrt(I3); //this is J=det(F)
    //derivatives of the reduced invariants 
    Eigen::Vector3d dIb = -1./3.*pow(I3, -4./3.)*I1*dI3 + pow(I3, -1./3.)*dI1;
    Eigen::Vector3d dJ = .5*pow(I3, -.5)*dI3;
    //second derivatives of the reduced invariants
    Eigen::Matrix3d d2Ib =  4./9.*pow(I3, -7./3.) * I1 * dI3*dI3.transpose()
                           -1./3.*pow(I3, -4./3.) * (dI1*dI3.transpose() + dI3*dI1.transpose())
                           -1./3.*pow(I3, -4./3.) * I1 * d2I3;
    Eigen::Matrix3d d2J =  -1./4.*pow(I3, -3./2.) * dI3*dI3.transpose() + .5/sqrt(I3)*d2I3;
    //obtaining the stresss and stiffness 
    Eigen::Vector3d s = _mu*dIb + 2*_K*(J-1)*dJ;
    result.stress << s(0), s(2), s(2), s(1);
    result.stiffness = 2*_mu*d2Ib +4*_K*dJ*dJ.transpose()  +4*_K*(J-1)*d2J;

  return stiffness;
  }
  
} //namespace TC
