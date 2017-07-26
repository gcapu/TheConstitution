#include "gtest/gtest.h"

TEST(otherdummy, suma)
{
    ASSERT_EQ(1+32, 33);
}

TEST(otherdummy, mult)
{
	int i = 11*3;
    ASSERT_EQ(i, 33);
}
