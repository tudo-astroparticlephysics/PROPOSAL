#include "gtest/gtest.h"
#include "PROPOSAL/propagation_utility/Decay.h"
#include "PROPOSAL/propagation_utility/DecayBuilder.h"
#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/crosssection/ParticleDefaultCrossSectionList.h"
//#include <memory>

using namespace PROPOSAL;

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}

auto GetCrossSections(bool interpol) {
    auto cuts = std::make_shared<EnergyCutSettings>(INF, 1, false);
    auto cross = GetStdCrossSections(MuMinusDef(), Ice(), cuts, interpol);
    return cross;
}

TEST(Constructor, Integral) {
    auto cross = GetCrossSections(false);
    auto decay = make_decay(cross, MuMinusDef(), false);
}

TEST(EnergyDecay, CompareIntegralInterpolant) {
    auto cross_interpol = GetCrossSections(true);
    auto decay_integral = make_decay(cross_interpol, MuMinusDef(), false);
    auto decay_interpol = make_decay(cross_interpol, MuMinusDef(), true);

    for (double logE_i = 3.; logE_i < 14; logE_i+=1) {
        double E_i = std::pow(10., logE_i);
        for (double rnd = 0; rnd < 1; rnd+=1e-1) {
            if (rnd >= 0.95)
                continue;
            double E_f_integral = decay_integral->EnergyDecay(E_i, rnd, Ice().GetMassDensity());
            double E_f_interpol = decay_interpol->EnergyDecay(E_i, rnd, Ice().GetMassDensity());
            EXPECT_NEAR(E_f_integral, E_f_interpol, E_f_integral*1e-3);
            EXPECT_GE(E_f_integral, MMU);
            EXPECT_GE(E_f_interpol, MMU);
            EXPECT_LE(E_f_integral, E_i);
            EXPECT_LE(E_f_interpol, E_i);
        }
        // test random numbers at the end of the phase space
        for (double rnd = 0.9999; rnd < 1; rnd+=1e-5) {
            if (rnd > 0.95)
                continue;
            double E_f_integral = decay_integral->EnergyDecay(E_i, rnd, Ice().GetMassDensity());
            double E_f_interpol = decay_interpol->EnergyDecay(E_i, rnd, Ice().GetMassDensity());
            EXPECT_NEAR(E_f_integral, E_f_interpol, E_f_integral*1e-3);
            EXPECT_GE(E_f_integral, MMU);
            EXPECT_GE(E_f_interpol, MMU);
            EXPECT_LE(E_f_integral, E_i);
            EXPECT_LE(E_f_interpol, E_i);
        }
    }

}