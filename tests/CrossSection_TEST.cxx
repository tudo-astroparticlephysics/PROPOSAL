#include "gtest/gtest.h"

#include "PROPOSAL/crosssection/parametrization/EpairProduction.h"
#include "PROPOSAL/crosssection/CrossSectionBuilder.h"

using namespace PROPOSAL;

TEST(CalculateStochasticLoss, GetUpperLimit_Exception)
{
    // This causes cubic_splines::find_parameter to fail in
    // CrossSectionDNDXInterpolant::GetUpperLimit, so we use Bisection as a
    // backup method
    auto param = crosssection::EpairForElectronPositron();
    auto medium = Air();
    auto cross = make_crosssection(
            param, EMinusDef(), medium,
            std::make_shared<EnergyCutSettings>(2.255, 1, false), true);

    double energy = 3443;
    double rate_failed = 0.000636779;
    auto comp_hash = medium.GetComponents().at(0).GetHash();

    // This function call causes the failure
    auto v_sampled = cross->CalculateStochasticLoss(comp_hash, energy,
                                                    rate_failed);

    // Check that value calculated by Bisection method is correct / consistent
    EXPECT_NEAR(cross->CalculateCumulativeCrosssection(energy, comp_hash, v_sampled),
                rate_failed, rate_failed*1e-5);
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
