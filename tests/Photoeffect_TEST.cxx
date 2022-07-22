#include "gtest/gtest.h"
#include "PROPOSAL/crosssection/parametrization/Photoeffect.h"
#include "PROPOSAL/crosssection/CrossSectionBuilder.h"
#include "PROPOSAL/medium/Medium.h"

using namespace PROPOSAL;

TEST(PhotoeffectSauter, dNdx_at_threshold)
{
    auto particle = GammaDef();
    auto medium = Air();
    auto photoeffect = crosssection::PhotoeffectSauter();

    auto cross = make_crosssection(photoeffect, particle, medium, nullptr, false);

    auto energy_lim = cross->GetLowerEnergyLim();

    // cross section should become zero below the energy_lim
    EXPECT_EQ(cross->CalculatedNdx(0.9*energy_lim), 0);
    // cross section should be zero at the energy_lim
    EXPECT_EQ(cross->CalculatedNdx(energy_lim), 0);
    // cross section should be non-zero above the energy_lim
    EXPECT_GT(cross->CalculatedNdx(1.001*energy_lim), 0);
}


TEST(PhotoeffectSauter, dEdx)
{
    auto particle = GammaDef();
    auto medium = Air();
    auto photoeffect = crosssection::PhotoeffectSauter();

    auto cross = make_crosssection(photoeffect, particle, medium, nullptr, false);

    EXPECT_EQ(cross->CalculatedEdx(1e3), 0);
    EXPECT_EQ(cross->CalculatedE2dx(1e3), 0);
}

TEST(PhotoEffectStauer, TestOnlyStochastic)
{
    auto particle = GammaDef();
    auto medium = Air();
    auto photoeffect = crosssection::PhotoeffectSauter();
    auto cuts = std::make_shared<EnergyCutSettings>(2, 1, false);

    // this should be an invalid configuration
    EXPECT_THROW(make_crosssection(photoeffect, particle, medium, cuts, false), std::logic_error);
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
