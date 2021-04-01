
#include "gtest/gtest.h"

#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/propagation_utility/PropagationUtility.h"
#include "PROPOSAL/propagation_utility/DisplacementBuilder.h"
#include "PROPOSAL/propagation_utility/TimeBuilder.h"
#include "PROPOSAL/propagation_utility/InteractionBuilder.h"
#include "PROPOSAL/crosssection/ParticleDefaultCrossSectionList.h"

using namespace PROPOSAL;

TEST(Comparison, Comparison_equal) {
    auto medium = Ice();
    auto cuts = std::make_shared<EnergyCutSettings>(500, 0.05, false);
    auto p_def = MuMinusDef();

    auto cross = GetStdCrossSections(p_def, medium, cuts, false);

    auto collection1 = PropagationUtility::Collection();
    collection1.displacement_calc = make_displacement(cross, false);
    collection1.interaction_calc = make_interaction(cross, false);
    collection1.time_calc = make_time(cross, p_def, false);

    auto collection2 = collection1;

    EXPECT_TRUE(collection1 == collection2);
}

TEST(Comparison, Comparison_not_equal) {
    auto medium = Ice();
    auto cuts = std::make_shared<EnergyCutSettings>(500, 0.05, false);
    auto p_def = MuMinusDef();

    auto cross = GetStdCrossSections(p_def, medium, cuts, false);

    auto collection1 = PropagationUtility::Collection();
    collection1.displacement_calc = make_displacement(cross, false);
    collection1.interaction_calc = make_interaction(cross, false);
    collection1.time_calc = make_time(cross, p_def, false);

    auto collection2 = PropagationUtility::Collection();
    collection2.displacement_calc = make_displacement(cross, false);
    collection2.interaction_calc = make_interaction(cross, false);
    collection2.time_calc = make_time(cross, p_def, false);

    EXPECT_FALSE(collection1 == collection2);
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
