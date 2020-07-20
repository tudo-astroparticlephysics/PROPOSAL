
#include "gtest/gtest.h"

#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/propagation_utility/PropagationUtility.h"

using namespace PROPOSAL;

TEST(Comparison, Comparison_equal) {
    auto water = std::make_shared<Water>(1.0);
    EnergyCutSettings ecuts;
    ParticleDef pDef(MuMinusDef::Get());
    PropagationUtility::Collection utility_defs;

    PropagationUtility utils1(pDef, water, ecuts, utility_defs);
    PropagationUtility utils2(pDef, water, ecuts, utility_defs);

    EXPECT_TRUE(utils1 == utils2);
}

TEST(Comparison, Comparison_not_equal) {
    auto water1 = std::make_shared<Water>(1.0);
    auto water2 = std::make_shared<Water>(0.9);

    EnergyCutSettings ecuts1(500, 0.05);
    EnergyCutSettings ecuts2(200, 0.01);

    ParticleDef pDef1(MuMinusDef::Get());
    ParticleDef pDef2(TauMinusDef::Get());

    PropagationUtility::Collection utility_defs;

    PropagationUtility utils1(pDef1, water1, ecuts1, utility_defs);
    PropagationUtility utils2(pDef2, water1, ecuts1, utility_defs);
    PropagationUtility utils3(pDef1, water2, ecuts1, utility_defs);
    PropagationUtility utils4(pDef1, water1, ecuts2, utility_defs);

    EXPECT_TRUE(utils1 != utils2);
    EXPECT_TRUE(utils1 != utils3);
    EXPECT_TRUE(utils1 != utils4);
}

TEST(Copyconstructor, Copyconstructor) {
    Utility A(MuMinusDef::Get(), std::make_shared<Ice>(), EnergyCutSettings(),
              Utility::Definition());
    Utility B(A);

    EXPECT_TRUE(A == B);

    Utility C(MuMinusDef::Get(), std::make_shared<Ice>(), EnergyCutSettings(),
              Utility::Definition(), InterpolationDef());
    Utility D(C);

    EXPECT_TRUE(C == D);
}

TEST(Copyconstructor, Copyconstructor2) {
    PropagationUtility A(MuMinusDef::Get(), std::make_shared<Medium>(Ice()), EnergyCutSettings(),
                         PropagationUtility::Collection());
    PropagationUtility B = A;

    EXPECT_TRUE(A == B);

    PropagationUtility C(MuMinusDef::Get(), std::make_shared<Medium>(Ice()), EnergyCutSettings(),
                         PropagationUtility::Collection(), InterpolationDef());
    PropagationUtility D = C;

    EXPECT_TRUE(C == D);
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
