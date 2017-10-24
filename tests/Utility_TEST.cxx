
#include "gtest/gtest.h"

#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/propagation_utility/PropagationUtility.h"

using namespace PROPOSAL;


TEST(Comparison , Comparison_equal )
{
    Water water1(1.0);
    Water water2(1.0);
    EXPECT_TRUE(water1 == water2);

    EnergyCutSettings ecuts1;
    EnergyCutSettings ecuts2;
    EXPECT_TRUE(ecuts1 == ecuts2);

    ParticleDef pDef1(MuMinusDef::Get());
    ParticleDef pDef2(MuMinusDef::Get());
    EXPECT_TRUE(pDef1 == pDef2);

    Utility::Definition utility_defs;

    Utility utils1(pDef1, water1, ecuts1, utility_defs);
    Utility utils2(pDef2, water2, ecuts2, utility_defs);

    EXPECT_TRUE(utils1==utils2);
}

int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
