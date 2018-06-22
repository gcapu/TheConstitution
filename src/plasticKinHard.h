#pragma once
#include <Eigen/Core>
#include "utils.h"
namespace TC
{

template<typename T> class traits;
template <typename _Scalar, int _Dim> class PlasticKinHard;

template<typename _Scalar, int _Dim> 
class traits<PlasticKinHard<_Scalar, _Dim>>{
 public:
  enum{
    Dim = _Dim,
    vDim = _Dim==3?6:4
  };
  using Scalar = _Scalar;
  using VectorType = Eigen::Matrix<_Scalar, vDim, 1>;
  using MatrixType = Eigen::Matrix<_Scalar, vDim, vDim>;
  class ResultType {
   public:
    ResultType(const VectorType& stress):_stress(stress){}
    const VectorType& stress() const {return _stress;}
   private:
    const VectorType& _stress;
  };
};

template <typename _Scalar, int _Dim> 
class PlasticKinHard{
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  using traits = traits<PlasticKinHard<_Scalar, _Dim>>;
  using Scalar = _Scalar;
  enum {
    Dim = _Dim,
    vDim = traits::vDim
  };
  using VectorType = typename traits::VectorType;
  using ResultType = typename traits::ResultType;
protected:
  //material properties
  const Scalar _lambda;
  const Scalar _mu;
  const Scalar _sigma_0;
  const Scalar _H;
  //internal variables
  Scalar plastic_strain;
  VectorType strain;
  VectorType stress;
  VectorType alpha;
public:
  PlasticKinHard(Scalar lambda, Scalar kappa, Scalar sigma_0, Scalar H);
  template<typename Derived>
  inline ResultType Increment(const Eigen::MatrixBase<Derived>& DStrain);
  //template<typename Derived>
  //inline ResultType Stiffness(const Eigen::MatrixBase<Derived>& F) const;
  };

template <typename _Scalar, int _Dim>
PlasticKinHard<_Scalar, _Dim>::PlasticKinHard(_Scalar lambda, _Scalar mu, _Scalar sigma_0, _Scalar hardening_slope): 
      _lambda(lambda), _mu(mu), _sigma_0(sigma_0), _H(hardening_slope), plastic_strain(0) {
  strain.setZero();
  stress.setZero();
  alpha.setZero();
  };

template <typename _Scalar, int _Dim>
template <typename Derived>
typename PlasticKinHard<_Scalar, _Dim>::ResultType
PlasticKinHard<_Scalar, _Dim>::Increment(const Eigen::MatrixBase<Derived>& DStrain) 
  {
  EIGEN_STATIC_ASSERT_SAME_MATRIX_SIZE(Derived, VectorType);
  VectorType _DStrain = DStrain.derived();
  Scalar dStrainTrace = vtrace(_DStrain); //magic constant is necessary since the universe is 3D
  // linear stress increment = lambda tr(DS) I + 2 mu DS
  stress += _lambda * dStrainTrace * videntity<Scalar, vDim>() + 2*_mu * _DStrain;
  //stress measured from the backstress (let's call it relative stress)
  VectorType relStress = stress - alpha; 
  //pressure p =
  Scalar p = vtrace(relStress)/3.;
  //deviatoric part of the relative stress 
  VectorType devRelStress = relStress - p*videntity<Scalar, vDim>(); 
  //the norm of the deviatoric part of the relative stress. Again it needs the out of plane component in 2D
  Scalar eta = vnorm(devRelStress);
  //obtaining the increment in gamma and the plastic strain
  Scalar diff = eta - sqrt(2./3.)*_sigma_0;
  Scalar cond = diff>0;
  Scalar dGamma = cond ? diff/(2*_mu+2./3.*_H) : 0; 
  plastic_strain += sqrt(2./3.)*dGamma;
  //updating the backstress
  Scalar factor = cond ? 2./3.*_H*dGamma/eta : 0;
  alpha+= factor * devRelStress;
  //updating the stress
  factor = cond ? 2.*_mu*dGamma/eta : 0;
  stress -= factor * devRelStress;
  return ResultType(stress);
  }
  
} //namespace TC
