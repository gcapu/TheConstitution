#include "gtest/gtest.h"
#include <Eigen/Core>
#include "plasticKinHard.h"
#include "utils.h"

#include <fstream>

template<typename DerivedA, typename DerivedB>
bool allclose(const Eigen::DenseBase<DerivedA>& a,
  const Eigen::DenseBase<DerivedB>& b,
  const typename DerivedA::RealScalar& rtol
  = Eigen::NumTraits<typename DerivedA::RealScalar>::dummy_precision(),
  const typename DerivedA::RealScalar& atol
  = Eigen::NumTraits<typename DerivedA::RealScalar>::epsilon())
  {
  return ((a.derived() - b.derived()).array().abs()
    <= (atol + rtol * b.derived().array().abs())).all();
  }

TEST(plasticKinematicHardening, doub3D_zero)
  {
  using matType = TC::PlasticKinHard<double, 3>;
  matType mat(2, 2, 1, 0.1);
  for(int i = 0; i<10; i++)
    mat.Increment(Eigen::Matrix<double,6,1>::Zero());
  matType::ResultType result = mat.Increment(Eigen::Matrix<double,6,1>::Zero());
  ASSERT_TRUE(allclose(result.stress(), Eigen::Matrix<double,6,1>::Zero()));
  }
TEST(plasticKinematicHardening, doub3D_zero_increment)
  {
  using matType = TC::PlasticKinHard<double, 3>;
  matType mat(2, 2, 1, 0.1);
  Eigen::Matrix<double, 6, 1> dstrain;
  dstrain << 0.1, 0, 0, 0, 0, 0, 0, 0, 0;
  for(int i = 0; i<10; i++)
    mat.Increment(dstrain);
  Eigen::Matrix<double, 6, 1> stress1 = mat.Increment(dstrain).stress();
  Eigen::Matrix<double, 6, 1> stress2 = mat.Increment(Eigen::Matrix<double,6,1>::Zero()).stress();
  ASSERT_GT(stress1(0,0), 0);
  ASSERT_TRUE(allclose(stress1, stress2));
  }

TEST(plasticKinematicHardening, doub3D_flows)
  {
  double E = 2, nu = 0;
  double lambda = TC::lambda(E, nu), mu = TC::mu(E, nu), sigma_0 = 1, H = 0;
  using matType = TC::PlasticKinHard<double, 3>;
  matType mat(lambda, mu, sigma_0, H);
  Eigen::Matrix<double,6,1> dstrain;
  double deVal=.1;
  dstrain << deVal, 0, 0, 0, 0, 0;
  for(int i = 0; i<10; i++)
    {
    Eigen::Matrix<double, 6, 1> stressLinear = mat.Increment(dstrain).stress();
    if(stressLinear(0) < sigma_0)
      ASSERT_DOUBLE_EQ(stressLinear(0), E*(i+1)*deVal);
    }
  Eigen::Matrix<double, 6, 1> stress1 = mat.Increment(dstrain).stress();
  Eigen::Matrix<double, 6, 1> dstress = mat.Increment(dstrain).stress()-stress1;
  ASSERT_LT(dstress(0,0), .5*(lambda+2*mu)*deVal);
  }

TEST(plasticKinematicHardening, float_3Dvs2D)
  {
  double E = 2, nu = 0.1;
  double lambda = TC::lambda(E, nu), mu = TC::mu(E, nu), sigma_0 = 1, H = 1;
  using matType3D = TC::PlasticKinHard<float, 3>;
  using matType2D = TC::PlasticKinHard<float, 2>;
  matType3D mat3D(lambda, mu, sigma_0, H);
  matType2D mat2D(lambda, mu, sigma_0, H);
  Eigen::Matrix<float,6,1> dstrain3D;
  Eigen::Matrix<float,4,1> dstrain2D;
  double deVal=.1;
  dstrain3D << deVal, 0, 0, 0, 0, 0;
  dstrain2D << deVal, 0, 0, 0;
  for(int i = 0; i<10; i++)
    {
    mat3D.Increment(dstrain3D);
    mat2D.Increment(dstrain2D);
    }
  Eigen::Matrix<float, 6, 1> stress3D = mat3D.Increment(dstrain3D).stress();
  Eigen::Matrix<float, 4, 1> stress2D = mat2D.Increment(dstrain2D).stress();
  ASSERT_FLOAT_EQ(stress3D(0), stress2D(0));
  ASSERT_FLOAT_EQ(stress3D(1), stress2D(1));
  ASSERT_FLOAT_EQ(stress3D(2), stress2D(2));
  ASSERT_FLOAT_EQ(stress3D(3), stress2D(3));
  for(int i = 0; i<25; i++)
    {
    mat3D.Increment(-dstrain3D);
    mat2D.Increment(-dstrain2D);
    }
  stress3D = mat3D.Increment(dstrain3D).stress();
  stress2D = mat2D.Increment(dstrain2D).stress();
  ASSERT_FLOAT_EQ(stress3D(0), stress2D(0));
  ASSERT_FLOAT_EQ(stress3D(1), stress2D(1));
  ASSERT_FLOAT_EQ(stress3D(2), stress2D(2));
  ASSERT_FLOAT_EQ(stress3D(3), stress2D(3));
  }
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
