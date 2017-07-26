#include "gtest/gtest.h"

TEST(subdummy, substraction)
{
    ASSERT_EQ(1, 2-1);
}

TEST(subdummy, mult)
{
	int i = 11*3;
    ASSERT_EQ(i, 33);
}
