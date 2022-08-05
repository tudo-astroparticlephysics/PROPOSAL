#include "gtest/gtest.h"
#include "PROPOSAL/crosssection/parametrization/Photoeffect.h"
#include "PROPOSAL/crosssection/parametrization/PhotoPairProduction.h"
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

    // should be zero for only-stochastic interaction
    EXPECT_EQ(cross->CalculatedEdx(1e3), 0);
    EXPECT_EQ(cross->CalculatedE2dx(1e3), 0);
}

TEST(PhotoEffectSauter, TestOnlyStochastic)
{
    auto particle = GammaDef();
    auto medium = Air();
    auto photoeffect = crosssection::PhotoeffectSauter();
    auto cuts = std::make_shared<EnergyCutSettings>(2, 1, false);

    // this should be an invalid configuration
    EXPECT_THROW(make_crosssection(photoeffect, particle, medium, cuts, false), std::logic_error);
}

TEST(PhotoEffectSauter, Comparison)
{
    auto particle = GammaDef();
    auto medium = Air();
    auto photoeffect = crosssection::PhotoeffectSauter();
    auto photopair = crosssection::PhotoPairKochMotz();

    auto cross_photoeffect = make_crosssection(photoeffect, particle, medium, nullptr, false);
    auto cross_photopair = make_crosssection(photopair, particle, medium, nullptr, false);

    // photoeffect should be dominant at low energies, but subdominant at high energies
    double low_E = 1e0;
    EXPECT_GT(cross_photoeffect->CalculatedNdx(low_E), cross_photopair->CalculatedNdx(low_E));

    double high_E = 1e6;
    EXPECT_LT(cross_photoeffect->CalculatedNdx(high_E), cross_photopair->CalculatedNdx(high_E));
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
