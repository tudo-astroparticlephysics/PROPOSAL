
#include "gtest/gtest.h"

#include <fstream>
#include "PROPOSAL/math/RandomGenerator.h"
#include "PROPOSAL/medium/Medium.h"
#include "PROPOSAL/medium/MediumFactory.h"
#include "PROPOSAL/particle/ParticleDef.h"
#include "PROPOSAL/propagation_utility/PropagationUtility.h"

using namespace PROPOSAL;

ParticleDef getParticleDef(const std::string& name) {
    if (name == "MuMinus") {
        return MuMinusDef::Get();
    } else if (name == "TauMinus") {
        return TauMinusDef::Get();
    } else {
        return EMinusDef::Get();
    }
}


TEST(ContinuousRandomization, Constructor) {
    auto cross_dummy = std::make_shared<CrossSectionBuilder>("contrand_dummy");
    cross_dummy->SetdEdx_function([](double energy)->double {return 5 + 1e-5 * energy;});
    cross_dummy->SetdE2dx_function([](double energy)->double {return std::exp(-14 + 2 * std::log(energy));});


    ContRand* dummy1 = new ContRandBuilder<UtilityIntegral>(CrossSectionList{cross_dummy}, getParticleDef("MuMinus"));
    ContRand* dummy2 = new ContRandBuilder<UtilityInterpolant>(CrossSectionList{cross_dummy}, getParticleDef("MuMinus"));

    auto dummy3 = ContRandBuilder<UtilityIntegral>(CrossSectionList{cross_dummy}, getParticleDef("MuMinus"));
    auto dummy4 = ContRandBuilder<UtilityInterpolant>(CrossSectionList{cross_dummy}, getParticleDef("MuMinus"));

}

class ContRandDummy : public ContRand{
public:
    ContRandDummy(CrossSectionList cross, const ParticleDef &def) : ContRand(cross, def){};
    double GetMass(){return mass;};
    double EnergyRandomize(double initial_energy, double final_energy, double rnd){ (void) initial_energy; (void) rnd; return final_energy;}
};

TEST(ContinuousRandomization, mass){
    auto cross_dummy = std::make_shared<CrossSectionBuilder>("contrand_dummy");
    cross_dummy->SetdEdx_function([](double energy)->double {return 5 + 1e-5 * energy;});
    cross_dummy->SetdE2dx_function([](double energy)->double {return std::exp(-14 + 2 * std::log(energy));});

    auto dummy1 = ContRandDummy(CrossSectionList{cross_dummy}, getParticleDef("MuMinus"));
    auto dummy2 = ContRandDummy(CrossSectionList{cross_dummy}, getParticleDef("TauMinus"));

    ASSERT_DOUBLE_EQ(dummy1.GetMass(), MMU);
    ASSERT_DOUBLE_EQ(dummy2.GetMass(), MTAU);
}

TEST(ContinuousRandomization, FirstMomentum){
    auto cross_dummy = std::make_shared<CrossSectionBuilder>("contrand_dummy");
    cross_dummy->SetdEdx_function([](double energy)->double {return 5 + 1e-5 * energy;});
    cross_dummy->SetdE2dx_function([](double energy)->double {return std::exp(-14 + 2 * std::log(energy));});

    RandomGenerator::Get().SetSeed(24601);
    auto contrand = ContRandBuilder<UtilityIntegral>(CrossSectionList{cross_dummy}, getParticleDef("MuMinus"));
    auto energies = std::array<std::pair<double, double>, 5>{{{5e4, 1e5}, {5e6, 1e7}, {5e8, 1e9}, {5e10, 1e11}, {5e12, 1e13}}};
    int statistics = 1e4;
    double average, randomized;
    for(auto E : energies){
        average = 0;
        for(int n=1; n<=statistics; n++){
            randomized = contrand.EnergyRandomize(E.second, E.first, RandomGenerator::Get().RandomDouble());
            average = average + (randomized - average) / n;
        }
        EXPECT_FALSE(average == E.first);
        EXPECT_NEAR(average, E.first, E.first * 1e-2);
    }
}

TEST(ContinuousRandomization, IdenticalEnergies){
    auto cross_dummy = std::make_shared<CrossSectionBuilder>("contrand_dummy");
    cross_dummy->SetdEdx_function([](double energy)->double {return 5 + 1e-5 * energy;});
    cross_dummy->SetdE2dx_function([](double energy)->double {return std::exp(-14 + 2 * std::log(energy));});
    auto contrand = ContRandBuilder<UtilityIntegral>(CrossSectionList{cross_dummy}, getParticleDef("MuMinus"));
    double energy = 1e5;
    double energy_randomized = contrand.EnergyRandomize(energy, energy, RandomGenerator::Get().RandomDouble());
    EXPECT_DOUBLE_EQ(energy_randomized, energy);
}

/*
TEST(ContinuousRandomization, Randomize_interpol) {
    std::string filename = "bin/TestFiles/continous_randomization.txt";
	std::ifstream in{filename};
	EXPECT_TRUE(in.good()) << "Test resource file '" << filename << "' could not be opened";


    char firstLine[256];
    in.getline(firstLine, 256);
    double initial_energy;
    double final_energy;
    double randomized_energy;
    double randomized_energy_new;
    double vcut;
    double ecut;
    std::string mediumName;
    std::string particleName;
    double rnd;
    double energy_old;
    bool first = true;

    std::cout.precision(16);

    RandomGenerator::Get().SetSeed(0);

    while (in.good()) {
        if (first)
            in >> rnd >> particleName >> mediumName >> ecut >> vcut >>
                initial_energy >> final_energy >> randomized_energy;

        first = false;
        energy_old = -1;

        std::shared_ptr<const Medium> medium           = CreateMedium(mediumName);
        EnergyCutSettings cut_settings(ecut, vcut);
        ParticleDef particle_def = getParticleDef(particleName);

        Utility utility(particle_def, medium, cut_settings,
                        Utility::Definition(), InterpolationDef());
        ContinuousRandomizer cont(utility, InterpolationDef());

        while (energy_old < initial_energy) {
            energy_old = initial_energy;

            randomized_energy_new =
                cont.Randomize(initial_energy, final_energy, rnd);

            ASSERT_NEAR(randomized_energy_new, randomized_energy,
                        1e-1 * randomized_energy);

            in >> rnd >> particleName >> mediumName >> ecut >> vcut >>
                initial_energy >> final_energy >> randomized_energy;
        }
    }
}
*/

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
