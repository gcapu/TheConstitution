#include "gtest/gtest.h"
#include "materialPoint.h"

TEST(dummy, cols)
  {
  Eigen::Matrix2f mat;
  ASSERT_EQ(cols(mat), 2);
  }

TEST(dummy, matrix)
  {
  Eigen::Matrix2f mat1 = Eigen::Matrix2f::Random();
  Eigen::Matrix2f mat2 = Eigen::Matrix2f::Random();
  Eigen::Matrix2f prod = mat1*mat2;
  Eigen::Matrix2f prod2 = matrixstuff(mat1,mat2);
  ASSERT_EQ(matrixstuff(mat1,mat2), prod);
  }

