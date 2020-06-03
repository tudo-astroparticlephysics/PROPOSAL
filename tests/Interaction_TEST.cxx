#include "gtest/gtest.h"

#include "PROPOSAL/propagation_utility/PropagationUtilityIntegral.h"
#include "PROPOSAL/propagation_utility/Interaction.h"
#include "PROPOSAL/crossection/CrossSectionBuilder.h"
#include "PROPOSAL/propagation_utility/PropagationUtility.h"
#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/math/RandomGenerator.h"


using namespace PROPOSAL;

TEST(Interaction, Constructor){
    auto cross = std::make_shared<CrossSectionBuilder>("interaction_dummy", MuMinusDef());
    cross->SetdEdx_function([](double energy)->double {return 5 + 1e-5 * energy;});
    cross->SetdNdx_function([](double energy)->double {
        //dNdx fit valid approximately for energies between 1e6 and 1e14 MeV
        double m = 0.03;
        double b = -11.;
        return std::exp(b + m * std::log(energy));
    });

    Interaction* interaction1 = new InteractionBuilder<UtilityIntegral>(CrossSectionList{cross});
    Interaction* interaction2 = new InteractionBuilder<UtilityInterpolant>(CrossSectionList{cross});

    InteractionBuilder<UtilityIntegral> interaction_3(CrossSectionList{cross});
    InteractionBuilder<UtilityInterpolant> interaction_4(CrossSectionList{cross});
}

class HistogramInteraction{
public:
    HistogramInteraction(){
        for(auto type : types){
            histogram[static_cast<int>(type)] = 0;
        }
        sum = 0;
    }
    void AddEntry(int type){ histogram[type] += 1; sum += 1;}
    InteractionType HighestCounter(){
        InteractionType max;
        int argmax = -1;
        for(auto type : types){
            if(histogram[static_cast<int>(type)] > argmax){
                max = type;
                argmax = histogram[static_cast<int>(type)];
            }
        }
        return max;
    }
    InteractionType LowestCounter(){
        InteractionType min;
        int argmin = std::numeric_limits<int>::max();
        for(auto type : types){
            if(histogram[static_cast<int>(type)] < argmin){
                min = type;
                argmin = histogram[static_cast<int>(type)];
            }
        }
        return min;
    }
private:
    const std::array<InteractionType, 4> types = {InteractionType::Brems, //Brems
                                      InteractionType::Ioniz, //Ioniz
                                      InteractionType::Epair, //Epair
                                      InteractionType::Photonuclear  //Photonuclear
                                      };
    std::unordered_map<int, int> histogram;
    double sum;
};

std::array<double, 2> GetRandomNumbers(){ return {RandomGenerator::Get().RandomDouble(), RandomGenerator::Get().RandomDouble()};}

TEST(TypeInteraction, Ratios){
    RandomGenerator::Get().SetSeed(24601);
    auto medium = std::make_shared<StandardRock>();
    auto ecuts = std::make_shared<EnergyCutSettings>(std::numeric_limits<double>::infinity(), 0.05, false);
    auto muon = MuMinusDef();
    auto cross = muon.GetCrossSections(medium, ecuts);

    auto interaction = InteractionBuilder<UtilityIntegral>(cross);
    HistogramInteraction histogram_low;
    HistogramInteraction histogram_high;
    int statistics = 1e3;
    for(int n=1; n<=statistics; n++){
        auto interaction_low = interaction.TypeInteraction(1e3, GetRandomNumbers());
        histogram_low.AddEntry(interaction_low->GetTypeId());

        auto interaction_high = interaction.TypeInteraction(1e10, GetRandomNumbers());
        histogram_high.AddEntry(interaction_high->GetTypeId());
    }

    EXPECT_EQ(histogram_low.HighestCounter(), InteractionType::Ioniz);
    EXPECT_EQ(histogram_high.LowestCounter(), InteractionType::Ioniz);

}

TEST(EnergyInteraction, Constraints){
    RandomGenerator::Get().SetSeed(24601);
    auto medium = std::make_shared<StandardRock>();
    auto ecuts = std::make_shared<EnergyCutSettings>(std::numeric_limits<double>::infinity(), 0.05, false);
    auto muon = MuMinusDef();
    auto cross = muon.GetCrossSections(medium, ecuts);

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

TEST(EnergyInteraction, compare_integral_interpolant1){
    RandomGenerator::Get().SetSeed(24601);
    auto medium = std::make_shared<StandardRock>();
    auto ecuts = std::make_shared<EnergyCutSettings>(std::numeric_limits<double>::infinity(), 0.05, false);
    auto muon = MuMinusDef();
    auto cross = muon.GetCrossSections(medium, ecuts);

    auto interaction_integral = InteractionBuilder<UtilityIntegral>(cross);
    auto interaction_interpol = InteractionBuilder<UtilityInterpolant>(cross);

    auto energies = std::array<double, 6>{1e3, 1e5, 1e7, 1e9, 1e11, 1e13};
    double rnd, energy_integral, energy_interpol;
    for(auto E_i : energies){
        rnd = RandomGenerator::Get().RandomDouble();
        energy_integral = interaction_integral.EnergyInteraction(E_i, rnd);
        energy_interpol = interaction_interpol.EnergyInteraction(E_i, rnd);
        EXPECT_NEAR(energy_integral, energy_interpol, energy_integral*1e-4);
    }
}

TEST(EnergyInteraction, compare_integral_interpolant2){
    RandomGenerator::Get().SetSeed(24601);
    auto medium = std::make_shared<StandardRock>();
    auto ecuts = std::make_shared<EnergyCutSettings>(500, 0.05, false);
    auto muon = MuMinusDef();
    auto cross = muon.GetCrossSections(medium, ecuts);

    auto interaction_integral = InteractionBuilder<UtilityIntegral>(cross);
    auto interaction_interpol = InteractionBuilder<UtilityInterpolant>(cross);

    auto energies = std::array<double, 6>{1e3, 1e5, 1e7, 1e9, 1e11, 1e13};
    double rnd, energy_integral, energy_interpol;
    for(auto E_i : energies){
        rnd = RandomGenerator::Get().RandomDouble();
        energy_integral = interaction_integral.EnergyInteraction(E_i, rnd);
        energy_interpol = interaction_interpol.EnergyInteraction(E_i, rnd);
        EXPECT_NEAR(energy_integral, energy_interpol, energy_integral*1e-4);
    }
}

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}