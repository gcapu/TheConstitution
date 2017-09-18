#pragma once
#include <Eigen/Core>

namespace TC
{

template <typename T>
class Material3D{
public:
  virtual Eigen::Matrix<T, 3, 3> Stress(const Eigen::Ref<const Eigen::Matrix<T, 3, 3>>& F) = 0;
  virtual Eigen::Matrix<T, 6, 6> Stiffness(const Eigen::Ref<const Eigen::Matrix<T, 3, 3>>& F) = 0;
  virtual ~Material3D(){}
  };

template <typename T>
class Material2D {
public:
  virtual Eigen::Matrix<T, 2, 2> Stress(const Eigen::Ref<const Eigen::Matrix<T, 2, 2>>& F) = 0;
  virtual Eigen::Matrix<T, 3, 3> Stiffness(const Eigen::Ref<const Eigen::Matrix<T, 2, 2>>& F) = 0;
  virtual ~Material2D() {}
  };

template <typename T>
class Material{
public:  
  virtual Eigen::Matrix<T, 2, 2> Stress(const Eigen::Ref<const Eigen::Matrix<T, 2, 2>>& F) = 0;
  virtual Eigen::Matrix<T, 3, 3> Stress(const Eigen::Ref<const Eigen::Matrix<T, 3, 3>>& F) = 0;
  virtual Eigen::Matrix<T, 3, 3> Stiffness(const Eigen::Ref<const Eigen::Matrix<T, 2, 2>>& F) = 0;
  virtual Eigen::Matrix<T, 6, 6> Stiffness(const Eigen::Ref<const Eigen::Matrix<T, 3, 3>>& F) = 0;
  virtual ~Material() {}
};
}//namespace TC
