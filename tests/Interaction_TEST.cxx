#include "gtest/gtest.h"

#include "PROPOSAL/Constants.h"
#include "PROPOSAL/propagation_utility/PropagationUtilityIntegral.h"
#include "PROPOSAL/propagation_utility/InteractionBuilder.h"
#include "PROPOSAL/crosssection/CrossSectionBuilder.h"
#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/math/RandomGenerator.h"
#include "PROPOSAL/crosssection/ParticleDefaultCrossSectionList.h"

using namespace PROPOSAL;

auto GetCrossSections() {
    auto cuts = std::make_shared<EnergyCutSettings>(INF, 0.05, false);
    static auto cross = GetStdCrossSections(MuMinusDef(), Ice(), cuts, true);
    return cross;
}

TEST(Interaction, Constructor){
    auto cross = GetCrossSections();
    Interaction* interaction1 = new InteractionBuilder<UtilityIntegral>(cross);
    Interaction* interaction2 = new InteractionBuilder<UtilityInterpolant>(cross);

    InteractionBuilder<UtilityIntegral> interaction_3(cross);
    InteractionBuilder<UtilityInterpolant> interaction_4(cross);

    delete interaction1;
    delete interaction2;
}

class HistogramInteraction{
public:
    HistogramInteraction(){
        for(auto type : types){
            histogram[type] = 0;
        }
        sum = 0;
    }
    void AddEntry(InteractionType type){ histogram[type] += 1; sum += 1;}
    InteractionType HighestCounter(){
        InteractionType max;
        int argmax = -1;
        for(auto type : types){
            if(histogram[type] > argmax){
                max = type;
                argmax = histogram[type];
            }
        }
        return max;
    }
    InteractionType LowestCounter(){
        InteractionType min;
        int argmin = std::numeric_limits<int>::max();
        for(auto type : types){
            if(histogram[type] < argmin){
                min = type;
                argmin = histogram[type];
            }
        }
        return min;
    }
private:
    const std::array<InteractionType, 4> types = {InteractionType::Brems,
                                                  InteractionType::Ioniz,
                                                  InteractionType::Epair,
                                                  InteractionType::Photonuclear
                                                 };
    std::unordered_map<InteractionType, int> histogram;
    double sum;
};

double rnd_number(){ return RandomGenerator::Get().RandomDouble();}

TEST(TypeInteraction, Ratios){
    // Ionization should be the dominant interaction type for low energies
    // and the most suppressed interaction type for high energies
    RandomGenerator::Get().SetSeed(24601);
    auto cross = GetCrossSections();
    auto interaction = InteractionBuilder<UtilityIntegral>(cross);
    HistogramInteraction histogram_low;
    HistogramInteraction histogram_high;
    int statistics = 1e3;

    for(int n=1; n<=statistics; n++){
        auto rates_low = interaction.Rates(1e3);
        auto interaction_low = interaction.SampleLoss(1e3, rates_low, rnd_number());
        histogram_low.AddEntry(std::get<0>(interaction_low));

        auto rates_high = interaction.Rates(1e10);
        auto interaction_high = interaction.SampleLoss(1e10, rates_high, rnd_number());
        histogram_high.AddEntry(std::get<0>(interaction_high));
    }

    EXPECT_EQ(histogram_low.HighestCounter(), InteractionType::Ioniz);
    EXPECT_EQ(histogram_high.LowestCounter(), InteractionType::Ioniz);
}

TEST(EnergyInteraction, Constraints){
    // sampled interaction energies should never be below the rest mass
    // and never be above the initial energy
    RandomGenerator::Get().SetSeed(24601);
    auto cross = GetCrossSections();
    auto interaction = InteractionBuilder<UtilityIntegral>(cross);
    int statistics = 1e3;

    auto energies = std::array<double, 5>{106, 1e4, 1e6, 1e8, 1e10};
    for(int i=0; i<statistics; i++){
        for(auto E_i : energies){
            double tmp = interaction.EnergyInteraction(E_i, RandomGenerator::Get().RandomDouble());
            EXPECT_TRUE(tmp >= MMU);
            EXPECT_TRUE(tmp <= E_i);
        }
    }
}

TEST(EnergyInteraction, CompareIntegralInterpolant){
    // comparing intergral and interpolant values
    RandomGenerator::Get().SetSeed(24601);
    auto cross = GetCrossSections();
    auto interaction_integral = InteractionBuilder<UtilityIntegral>(cross);
    auto interaction_interpol = InteractionBuilder<UtilityInterpolant>(cross);

    double rnd, energy_integral, energy_interpol;
    for(double Elog_i = 3; Elog_i < 10; Elog_i+=1e-2 ){
        double E_i = std::pow(Elog_i, 10.);
        rnd = RandomGenerator::Get().RandomDouble();
        energy_integral = interaction_integral.EnergyInteraction(E_i, rnd);
        energy_interpol = interaction_interpol.EnergyInteraction(E_i, rnd);
        EXPECT_NEAR(energy_integral, energy_interpol, energy_integral * 1e-4);
    }
}

TEST(MeanFreePath, ConsistencyCheck)
{
    // The free mean path length should decrease for higher energies
    auto cross = GetCrossSections();
    auto interaction = InteractionBuilder<UtilityIntegral>(cross);

    double pathlength_old = INF;
    for(double Elog = 6; Elog < 10; Elog+=1e-3) {
        double E = std::pow(Elog, 10.);
        double pathlength = interaction.MeanFreePath(E);
        EXPECT_LT(pathlength, pathlength_old);
        pathlength_old = pathlength;
    }
}

TEST(MeanFreePath, CompareIntegralInterpolant)
{
    // comparing intergral and interpolant values
    auto cross = GetCrossSections();
    auto inter_integral = InteractionBuilder<UtilityIntegral>(cross);
    auto inter_interpol = InteractionBuilder<UtilityInterpolant>(cross);

    for(double Elog = 3; Elog < 10; Elog+=1e-3) {
        double E = std::pow(Elog, 10.);
        double pathlength_integral = inter_integral.MeanFreePath(E);
        double pathlength_interpol = inter_interpol.MeanFreePath(E);
        EXPECT_NEAR(pathlength_integral, pathlength_interpol, pathlength_integral * 1e-3);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
