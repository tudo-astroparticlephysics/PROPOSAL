
#include "gtest/gtest.h"

#include "PROPOSAL/particle/ParticleDef.h"

using namespace PROPOSAL;

TEST(Comparison, Comparison_equal)
{
    MuMinusDef A;
    MuMinusDef B;
    EXPECT_TRUE(A == B);

    ParticleDef* C = new MuMinusDef();
    ParticleDef* D = new MuMinusDef();

    EXPECT_TRUE(*C == *D);
    delete C;
    delete D;

    C = new ParticleDef(MuMinusDef());
    D = new ParticleDef(MuMinusDef());

    EXPECT_TRUE(*C == *D);
    delete C;
    delete D;

    C = new ParticleDef(TauMinusDef());
    D = new ParticleDef(TauMinusDef());

    EXPECT_TRUE(*C == *D);
    delete C;
    delete D;
}

TEST(Comparison, Comparison_not_equal)
{
    MuMinusDef A;
    MuMinusDef B;
    EXPECT_TRUE(A == B);

    ParticleDef AA = ParticleDef::Builder().SetMass(100).build();
    EXPECT_TRUE(AA != B);

    ParticleDef* C = new ParticleDef(MuMinusDef());
    ParticleDef* D = new ParticleDef(TauMinusDef());

    EXPECT_TRUE(*C != *D);
    delete C;
    delete D;
}

TEST(Assignment, Copyconstructor)
{
    MuMinusDef A;
    ParticleDef B = A;

    EXPECT_TRUE(A == B);
}

TEST(Assignment, Copyconstructor2)
{
    auto A = MuMinusDef();
    ParticleDef B(A);
    EXPECT_TRUE(A == B);
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
