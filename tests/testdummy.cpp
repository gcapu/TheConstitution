#include "gtest/gtest.h"
#include "materialBase.h"

TEST(dummy, cols)
  {
  Eigen::Matrix2f mat;
  ASSERT_EQ(TC::cols(mat), 2);
  }

