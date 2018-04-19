
#include "gtest/gtest.h"

#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/propagation_utility/PropagationUtility.h"

using namespace PROPOSAL;

TEST(Comparison, Comparison_equal)
{
    Water water(1.0);
    EnergyCutSettings ecuts;
    ParticleDef pDef(MuMinusDef::Get());
    Utility::Definition utility_defs;

    Utility utils1(pDef, water, ecuts, utility_defs);
    Utility utils2(pDef, water, ecuts, utility_defs);

    EXPECT_TRUE(utils1 == utils2);
}

TEST(Comparison, Comparison_not_equal)
{
    Water water1(1.0);
    Water water2(0.9);

    EnergyCutSettings ecuts1(500, 0.05);
    EnergyCutSettings ecuts2(200, 0.01);

    ParticleDef pDef1(MuMinusDef::Get());
    ParticleDef pDef2(TauMinusDef::Get());

    Utility::Definition utility_defs;

    Utility utils1(pDef1, water1, ecuts1, utility_defs);
    Utility utils2(pDef2, water1, ecuts1, utility_defs);
    Utility utils3(pDef1, water2, ecuts1, utility_defs);
    Utility utils4(pDef1, water1, ecuts2, utility_defs);

    EXPECT_TRUE(utils1 != utils2);
    EXPECT_TRUE(utils1 != utils3);
    EXPECT_TRUE(utils1 != utils4);
}

TEST(Copyconstructor, Copyconstructor)
{
    Utility A(MuMinusDef::Get(), Ice(), EnergyCutSettings(), Utility::Definition());
    Utility B(A);

    EXPECT_TRUE(A == B);

    Utility C(MuMinusDef::Get(), Ice(), EnergyCutSettings(), Utility::Definition(), InterpolationDef());
    Utility D(C);

    EXPECT_TRUE(C == D);
}

TEST(Copyconstructor, Copyconstructor2)
{
    Utility A(MuMinusDef::Get(), Ice(), EnergyCutSettings(), Utility::Definition());
    Utility B = A;

    EXPECT_TRUE(A == B);

    Utility C(MuMinusDef::Get(), Ice(), EnergyCutSettings(), Utility::Definition(), InterpolationDef());
    Utility D = C;

    EXPECT_TRUE(C == D);
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
