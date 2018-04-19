
#include "gtest/gtest.h"

#include "PROPOSAL/particle/ParticleDef.h"

using namespace PROPOSAL;

TEST(Comparison, Comparison_equal)
{
    ParticleDef A;
    ParticleDef B;
    EXPECT_TRUE(A == B);

    ParticleDef* C = new ParticleDef();
    ParticleDef* D = new ParticleDef();

    EXPECT_TRUE(*C == *D);
    delete C;
    delete D;

    C = new ParticleDef(MuMinusDef::Get());
    D = new ParticleDef(MuMinusDef::Get());

    EXPECT_TRUE(*C == *D);
    delete C;
    delete D;

    C = new ParticleDef(TauMinusDef::Get());
    D = new ParticleDef(TauMinusDef::Get());

    EXPECT_TRUE(*C == *D);
    delete C;
    delete D;
}

TEST(Comparison, Comparison_not_equal)
{
    ParticleDef A;
    ParticleDef B;
    EXPECT_TRUE(A == B);

    ParticleDef AA = ParticleDef::Builder().SetMass(100).build();
    EXPECT_TRUE(AA != B);

    ParticleDef* C = new ParticleDef(MuMinusDef::Get());
    ParticleDef* D = new ParticleDef(TauMinusDef::Get());

    EXPECT_TRUE(*C != *D);
    delete C;
    delete D;
}

TEST(Assignment, Copyconstructor)
{
    ParticleDef A;
    ParticleDef B = A;

    EXPECT_TRUE(A == B);
}

TEST(Assignment, Copyconstructor2)
{
    ParticleDef A(MuMinusDef::Get());
    ParticleDef B(A);
    EXPECT_TRUE(A == B);
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
