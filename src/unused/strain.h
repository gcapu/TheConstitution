#pragma once
#include <Eigen/Core>

namespace TC
{
template <typename _Scalar, int _Dim>
class StrainBase {
public:
  typedef _Scalar Scalar;
  enum {
    Dim = _Dim,
  typedef Eigen::Matrix<Scalar, Dim, Dim> MatrixType;
private:
  MatrixType _strain;
public:
  StrainBase(const MatrixType& strain) {_strain = strain;}
  operator const MatrixType&() {return _strain;}
  operator MatrixType() {return _strain;}
};

template <typename _Scalar, int _Dim, bool _symmetric = true>
class Strain;

template <typename _Scalar>
class Strain<_Scalar, 3, true>: StrainBase<_Scalar, 3>{
  typedef _Scalar Scalar;
  enum {
    Dim = 3,
    VoigtDim = 6
    };
  typedef Eigen::Matrix<Scalar, VoigtDim, 1> VoigtType;
  Strain(const MatrixType& strain): StrainBase(strain){}
  VoigtType ToVoigt(){
    VoigtType V;
    V << _strain(0,0), _strain(1,1), _strain(2,2), 
         2*_strain(0,1), 2*_strain(1,2), 2*_strain(0,2);
    return V;
  }
}

template <typename _Scalar>
class Strain<_Scalar, 3, false>: StrainBase<_Scalar, 3>{
  typedef _Scalar Scalar;
  enum {
    Dim = 3,
    VoigtDim = 9
    };
  typedef Eigen::Matrix<Scalar, VoigtDim, 1> VoigtType;
  Strain(const MatrixType& strain): StrainBase(strain){}
  VoigtType ToVoigt(){
    VoigtType V;
    V << _strain(0,0), _strain(1,1), _strain(2,2), 
         2*_strain(0,1), 2*_strain(1,2), 2*_strain(0,2),
         2*_strain(1,0), 2*_strain(2,1), 2*_strain(2,0);
    return V;
  }
}

template <typename _Scalar>
class Strain<_Scalar, 2, true>: StrainBase<_Scalar, 2>{
  typedef _Scalar Scalar;
  enum {
    Dim = 2,
    VoigtDim = 3
    };
  typedef Eigen::Matrix<Scalar, VoigtDim, 1> VoigtType;
  Strain(const MatrixType& strain): StrainBase(strain){}
  VoigtType ToVoigt(){
    VoigtType V;
    V << _strain(0,0), _strain(1,1), _strain(0,1);
    return V;
  }
}

template <typename _Scalar>
class Strain<_Scalar, 2, false>: StrainBase<_Scalar, 2>{
  typedef _Scalar Scalar;
  enum {
    Dim = 2,
    VoigtDim = 3
    };
  typedef Eigen::Matrix<Scalar, VoigtDim, 1> VoigtType;
  Strain(const MatrixType& strain): StrainBase(strain){}
  VoigtType ToVoigt(){
    VoigtType V;
    V << _strain(0,0), _strain(1,1), _strain(0,1), _strain(1,0);
    return V;
  }
}

}//namespace TC




