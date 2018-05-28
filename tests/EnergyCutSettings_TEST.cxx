
#include "gtest/gtest.h"

#include "PROPOSAL/EnergyCutSettings.h"

using namespace PROPOSAL;

TEST(Comparison, Comparison_equal)
{
    EnergyCutSettings A;
    EnergyCutSettings B;
    EXPECT_TRUE(A == B);
    EnergyCutSettings* C = new EnergyCutSettings(100, 0.01);
    EnergyCutSettings* D = new EnergyCutSettings(100, 0.01);
    EXPECT_TRUE(*C == *D);
    EnergyCutSettings* E = new EnergyCutSettings(500, 0.05);
    EXPECT_TRUE(A == *E);
}

TEST(Comparison, Comparison_not_equal)
{
    EnergyCutSettings A;
    EnergyCutSettings B(200, 0.09);
    EXPECT_TRUE(A != B);
    EnergyCutSettings* C = new EnergyCutSettings(200, 0.01);
    EnergyCutSettings* D = new EnergyCutSettings(100, 0.01);
    EXPECT_TRUE(*C != *D);
}
TEST(Assignment, Copyconstructor)
{
    EnergyCutSettings A;
    EnergyCutSettings B = A;

    EXPECT_TRUE(A == B);
}

TEST(Assignment, Copyconstructor2)
{
    EnergyCutSettings A(5000, 0.1);
    EnergyCutSettings B(A);

    EXPECT_TRUE(A == B);
}

TEST(Assignment, Operator)
{
    EnergyCutSettings A;
    EnergyCutSettings B(200, 0.01);

    EXPECT_TRUE(A != B);

    B = A;

    EXPECT_TRUE(A == B);

    A.SetEcut(300);

    EXPECT_TRUE(A != B);

    B = A;

    EXPECT_TRUE(A == B);
}

TEST(Assignment, Swap)
{
    EnergyCutSettings A;
    EnergyCutSettings B;
    EXPECT_TRUE(A == B);
    EnergyCutSettings* C = new EnergyCutSettings(100, 0.01);
    EnergyCutSettings* D = new EnergyCutSettings(100, 0.01);
    EXPECT_TRUE(*C == *D);
    A.swap(*C);
    EXPECT_TRUE(A == *D);
    EXPECT_TRUE(B == *C);
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
