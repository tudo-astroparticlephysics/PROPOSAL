#include "gtest/gtest.h"

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/crosssection/CrossSectionBuilder.h"
#include "PROPOSAL/crosssection/ParticleDefaultCrossSectionList.h"
#include "PROPOSAL/math/RandomGenerator.h"
#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/propagation_utility/InteractionBuilder.h"
#include "PROPOSAL/propagation_utility/PropagationUtilityIntegral.h"
#include "PROPOSAL/propagation_utility/DisplacementBuilder.h"

#include <memory>

using namespace PROPOSAL;

auto GetCrossSections()
{
    auto cuts = std::make_shared<EnergyCutSettings>(INF, 0.05, false);
    static auto cross = GetStdCrossSections(MuMinusDef(), Ice(), cuts, true);
    return cross;
}

TEST(Interaction, Constructor)
{
    auto cross = GetCrossSections();
    auto disp = std::shared_ptr<Displacement>(make_displacement(cross, false));
    auto interaction1
        = std::make_unique<InteractionBuilder>(disp, cross, std::false_type());
    auto interaction2
        = std::make_unique<InteractionBuilder>(disp, cross, std::true_type());

    InteractionBuilder interaction_3(disp, cross, std::false_type());
    InteractionBuilder interaction_4(disp, cross, std::true_type());
}

class HistogramInteraction {
public:
    HistogramInteraction()
    {
        for (auto type : types) {
            histogram[type] = 0;
        }
        sum = 0;
    }
    void AddEntry(InteractionType type)
    {
        histogram[type] += 1;
        sum += 1;
    }
    InteractionType HighestCounter()
    {
        InteractionType max;
        int argmax = -1;
        for (auto type : types) {
            if (histogram[type] > argmax) {
                max = type;
                argmax = histogram[type];
            }
        }
        return max;
    }
    InteractionType LowestCounter()
    {
        InteractionType min;
        int argmin = std::numeric_limits<int>::max();
        for (auto type : types) {
            if (histogram[type] < argmin) {
                min = type;
                argmin = histogram[type];
            }
        }
        return min;
    }

private:
    const std::array<InteractionType, 4> types
        = { InteractionType::Brems, InteractionType::Ioniz,
              InteractionType::Epair, InteractionType::Photonuclear };
    std::unordered_map<InteractionType, int> histogram;
    double sum;
};

double rnd_number() { return RandomGenerator::Get().RandomDouble(); }

TEST(TypeInteraction, Ratios)
{
    // Ionization should be the dominant interaction type for low energies
    // and the most suppressed interaction type for high energies
    RandomGenerator::Get().SetSeed(24601);
    auto cross = GetCrossSections();
    auto interaction = make_interaction(cross, false);
    HistogramInteraction histogram_low;
    HistogramInteraction histogram_high;
    int statistics = 1e3;

    for (int n = 1; n <= statistics; n++) {
        auto rates_low = interaction->Rates(1e3);
        auto interaction_low
            = interaction->SampleLoss(1e3, rates_low, rnd_number());
        histogram_low.AddEntry(interaction_low.type);

        auto rates_high = interaction->Rates(1e10);
        auto interaction_high
            = interaction->SampleLoss(1e10, rates_high, rnd_number());
        histogram_high.AddEntry(interaction_high.type);
    }

    EXPECT_EQ(histogram_low.HighestCounter(), InteractionType::Ioniz);
    EXPECT_EQ(histogram_high.LowestCounter(), InteractionType::Ioniz);
}

TEST(EnergyInteraction, Constraints)
{
    // sampled interaction energies should never be below the rest mass
    // and never be above the initial energy
    RandomGenerator::Get().SetSeed(24601);
    auto cross = GetCrossSections();
    auto interaction = make_interaction(cross, false);
    int statistics = 1e3;

    auto energies = std::array<double, 5> { 200, 1e4, 1e6, 1e8, 1e10 };
    for (int i = 0; i < statistics; i++) {
        for (auto E_i : energies) {
            double tmp = interaction->EnergyInteraction(
                E_i, RandomGenerator::Get().RandomDouble());
            EXPECT_TRUE(tmp >= MMU);
            EXPECT_TRUE(tmp <= E_i);
        }
    }
}

TEST(EnergyInteraction, CompareIntegralInterpolant)
{
    // comparing integral and interpolant values
    RandomGenerator::Get().SetSeed(24601);
    auto cross = GetCrossSections();
    auto low = CrossSectionVector::GetLowerLim(cross);
    auto interaction_integral = make_interaction(cross, false);
    auto interaction_interpol = make_interaction(cross, true);

    double rnd, energy_integral, energy_interpol;
    for (double Elog_i = std::log10(low); Elog_i < 14; Elog_i += 5e-2) {
        double E_i = std::pow(10., Elog_i);
        for (auto i_stat = 0; i_stat < 100; i_stat++) {
            rnd = RandomGenerator::Get().RandomDouble();
            energy_integral = interaction_integral->EnergyInteraction(E_i, rnd);
            energy_interpol = interaction_interpol->EnergyInteraction(E_i, rnd);
            auto lower_lim = CrossSectionVector::GetLowerLim(cross);
            auto rnd_integral = std::min(interaction_integral->EnergyIntegral(E_i, lower_lim), -std::log(rnd));
            auto rnd_interpol = std::min(interaction_interpol->EnergyIntegral(E_i, lower_lim), -std::log(rnd));
            double precision = 1e-3;
            if (E_i < 1e5)
                precision = 5e-3; // integrand hard to interpolate
            else if (E_i < 1e6)
                precision = 2e-3;
            EXPECT_NEAR(energy_integral, energy_interpol, energy_integral * precision);
            EXPECT_NEAR(interaction_integral->EnergyIntegral(E_i, energy_integral), rnd_integral, 1e-5);
            EXPECT_NEAR(interaction_interpol->EnergyIntegral(E_i, energy_interpol), rnd_interpol, 1e-5);
        }
    }
}

TEST(MeanFreePath, ConsistencyCheck)
{
    // The free mean path length should decrease for higher energies
    auto cross = GetCrossSections();
    auto interaction = make_interaction(cross, false);

    double pathlength_old = INF;
    for (double Elog = 6; Elog < 10; Elog += 1e-3) {
        double E = std::pow(Elog, 10.);
        double pathlength = interaction->MeanFreePath(E);
        EXPECT_LT(pathlength, pathlength_old);
        pathlength_old = pathlength;
    }
}

TEST(MeanFreePath, CompareIntegralInterpolant)
{
    // comparing intergral and interpolant values
    auto cross = GetCrossSections();
    auto inter_integral = make_interaction(cross, false);
    auto inter_interpol = make_interaction(cross, true);

    for (double Elog = 3; Elog < 10; Elog += 1e-3) {
        double E = std::pow(Elog, 10.);
        double pathlength_integral = inter_integral->MeanFreePath(E);
        double pathlength_interpol = inter_interpol->MeanFreePath(E);
        EXPECT_NEAR(pathlength_integral, pathlength_interpol,
            pathlength_integral * 1e-3);
    }
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
